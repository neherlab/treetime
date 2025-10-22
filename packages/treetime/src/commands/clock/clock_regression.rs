use super::clock_model::ClockModel;
use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::commands::clock::find_best_root::find_best_root::find_best_root;
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::reroot_in_place;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use clap::Args;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;

#[derive(Debug, Clone, Serialize, Deserialize, Args, SmartDefault)]
pub struct ClockOptions {
  /// Variance scaling factor proportional to branch length
  #[clap(long, default_value_t = ClockOptions::default().variance_factor)]
  #[default = 0.0]
  pub variance_factor: f64,

  /// Constant variance offset for all branches
  #[clap(long, default_value_t = ClockOptions::default().variance_offset)]
  #[default = 0.0]
  pub variance_offset: f64,

  /// Additional variance offset for leaf (terminal) nodes
  #[clap(long, default_value_t = ClockOptions::default().variance_offset_leaf)]
  #[default = 1.0]
  pub variance_offset_leaf: f64,
}

pub fn root_clock_model<N, E, D>(graph: &Graph<N, E, D>) -> Result<ClockModel, Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let root = graph.get_exactly_one_root()?;
  let root = root.read_arc().payload().read_arc();
  ClockModel::new(root.clock_set())
}

pub fn clock_regression_backward<N, E, D>(graph: &Graph<N, E, D>, options: &ClockOptions)
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  graph.par_iter_breadth_first_backward(|mut n| {
    let date = n.payload.date();
    let q_to_parent = if n.is_leaf {
      if n.payload.is_outlier() {
        ClockSet::outlier_contribution()
      } else {
        ClockSet::leaf_contribution(date)
      }
    } else {
      let mut q_dest = ClockSet::default();
      for (_, e) in &n.children {
        q_dest += e.read_arc().from_child();
      }
      q_dest
    };

    let is_leaf = n.is_leaf;
    let is_root = n.is_root;
    let edge_to_parent = n.get_exactly_one_parent_edge();

    if is_root {
      *n.payload.clock_set_mut() = q_to_parent;
    } else {
      let edge_to_parent = edge_to_parent.expect("Encountered a node without a parent edge");
      *edge_to_parent.to_parent_mut() = q_to_parent;
      let edge_len = edge_to_parent.branch_length().expect("Encountered an edge without a weight");
      let mut branch_variance = options.variance_factor * edge_len + options.variance_offset;
      *edge_to_parent.from_child_mut() = if is_leaf {
        branch_variance += options.variance_offset_leaf;
        ClockSet::leaf_contribution_to_parent(date, edge_len, branch_variance)
      } else {
        edge_to_parent.to_parent().propagate_averages(edge_len, branch_variance)
      };
    }
    GraphTraversalContinuation::Continue
  });
}

pub fn clock_regression_forward<N, E, D>(graph: &Graph<N, E, D>, options: &ClockOptions)
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Sync + Send,
{
  graph.par_iter_breadth_first_forward(|mut n| {
    if !n.is_root {
      let (_, edge) = n.get_exactly_one_parent().unwrap();
      let edge = edge.read_arc();
      let edge_len = edge.branch_length().expect("Encountered an edge without a weight");
      let branch_variance = options.variance_factor * edge_len + options.variance_offset;

      let mut q_dest = edge.to_parent().clone();
      q_dest += edge.to_child().propagate_averages(edge_len, branch_variance);
      *n.payload.clock_set_mut() = q_dest;
    }

    for mut child_edge in n.child_edges {
      let mut q = n.payload.clock_set().clone();
      q -= child_edge.from_child();
      *child_edge.to_child_mut() = q;
    }
    GraphTraversalContinuation::Continue
  });
}

pub fn estimate_clock_model_with_reroot<N, E, D>(
  graph: &mut Graph<N, E, D>,
  options: &ClockOptions,
  keep_root: bool,
  optimization_params: &BranchPointOptimizationParams,
) -> Result<ClockModel, Report>
where
  N: GraphNode + ClockNode + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  log::info!("## Estimating clock model (keep_root={})", keep_root);

  log::info!("### Running backward regression");
  clock_regression_backward(graph, options);
  log::debug!("Backward regression completed");

  if !keep_root {
    log::info!("### Running forward regression to find optimal root");
    clock_regression_forward(graph, options);
    log::debug!("Forward regression completed");

    log::info!("### Finding best root and rerooting tree");
    let new_root_key = reroot_in_place(graph, options, optimization_params)?;
    log::info!("Rerooted to node {}", new_root_key.0);
    log::debug!("Rerooting completed");
  } else {
    log::info!("### Keeping original root (--keep-root enabled)");
  }

  log::info!("### Extracting clock model from root");
  let clock_model = root_clock_model(graph)?;
  log::info!("**Clock rate:** {:.6e}", clock_model.clock_rate());
  log::info!("**Intercept:** {:.4}", clock_model.intercept());
  log::info!("**R²:** {:.4}", clock_model.r_val() * clock_model.r_val());
  log::debug!("Clock model: {}", serde_json::to_string_pretty(&clock_model).unwrap_or_else(|_| "<serialization failed>".to_owned()));

  Ok(clock_model)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::clock::clock_graph::ClockGraph;
  use crate::graph::node::Named;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::seq::div::{OnlyLeaves, calculate_divs};
  use approx::assert_ulps_eq;
  use maplit::btreemap;
  use std::collections::BTreeMap;

  pub fn calculate_naive_rate(dates: &BTreeMap<String, f64>, div: &BTreeMap<String, f64>) -> f64 {
    let t: f64 = dates.values().sum();
    let tsq: f64 = dates.values().map(|&x| x * x).sum();
    let dt: f64 = div.iter().map(|(c, div)| div * dates[c]).sum();
    let d: f64 = div.values().sum();
    (dt * 4.0 - d * t) / (tsq * 4.0 - (t * t))
  }

  #[test]
  fn test_clock_naive_rate() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let graph: ClockGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let divs = calculate_divs(&graph, OnlyLeaves(true));
    let naive_rate = calculate_naive_rate(&dates, &divs);

    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().date = Some(dates[&name]);
    }

    clock_regression_backward(&graph, &ClockOptions::default());
    let clock = root_clock_model(&graph)?;
    assert_ulps_eq!(naive_rate, clock.clock_rate(), epsilon = 1e-9);

    let options = &ClockOptions {
      variance_factor: 1.0,
      variance_offset: 0.0,
      variance_offset_leaf: 1.0,
    };

    clock_regression_backward(&graph, options);
    let clock = root_clock_model(&graph)?;
    assert_ulps_eq!(0.007710610998647367, clock.clock_rate(), epsilon = 1e-9);

    Ok(())
  }
}
