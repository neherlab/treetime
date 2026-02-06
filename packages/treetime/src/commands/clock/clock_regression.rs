use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_set::ClockSet;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::reroot_in_place;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use clap::Args;
use eyre::Report;
use log::{debug, info};
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

pub fn clock_regression_backward<N, E, D>(graph: &Graph<N, E, D>, options: &ClockOptions)
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  graph.par_iter_breadth_first_backward(|mut n| {
    let date = n.payload.likely_time();
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
      let edge_len = edge_to_parent
        .branch_length()
        .expect("Encountered an edge without a weight");
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

/// Estimates clock model with optional rerooting
pub fn estimate_clock_model_with_reroot<N, E, D>(
  graph: &mut Graph<N, E, D>,
  options: &ClockOptions,
  clock_rate: Option<f64>,
  keep_root: bool,
  optimization_params: &BranchPointOptimizationParams,
) -> Result<ClockModel, Report>
where
  N: GraphNode + ClockNode + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  if let Some(rate) = clock_rate {
    info!("## Estimating clock model with fixed rate {rate:.6e} (keep_root={keep_root})");
  } else {
    info!("## Estimating clock model (keep_root={keep_root})");
  }

  info!("### Running backward regression");
  clock_regression_backward(graph, options);
  debug!("Backward regression completed");

  if !keep_root {
    info!("### Running forward regression to find optimal root");
    clock_regression_forward(graph, options);
    debug!("Forward regression completed");

    info!("### Finding best root and rerooting tree");
    let new_root_key = reroot_in_place(graph, options, optimization_params)?;
    info!("Rerooted to node {}", new_root_key.0);
    debug!("Rerooting completed");
  } else {
    info!("### Keeping original root (--keep-root enabled)");
  }

  info!("### Extracting clock model from root");
  let root = graph.get_exactly_one_root()?;
  let root = root.read_arc().payload().read_arc();
  let clock_model = if let Some(rate) = clock_rate {
    info!("### Using fixed clock rate: {rate:.6e}");
    ClockModel::with_fixed_rate(root.clock_set(), rate)
  } else {
    info!("### Using estimated clock rate");
    ClockModel::new(root.clock_set())?
  };

  info!("**Clock rate:** {:.6e}", clock_model.clock_rate());
  info!("**Intercept:** {:.4}", clock_model.intercept());
  if let Some(r_val) = clock_model.r_val() {
    info!("**R²:** {:.4}", r_val * r_val);
  }
  if let Some(chisq) = clock_model.chisq() {
    info!("**χ²:** {chisq:.4}");
  }
  if let Some(hessian) = clock_model.hessian() {
    info!("**Hessian:**\n{hessian}");
  }
  debug!(
    "Clock model: {}",
    serde_json::to_string_pretty(&clock_model).unwrap_or_else(|_| "<serialization failed>".to_owned())
  );

  Ok(clock_model)
}
