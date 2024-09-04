use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_set::ClockSet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use super::clock_model::ClockModel;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClockOptions {
  pub variance_factor: f64,
  pub variance_offset: f64,
  pub variance_offset_leaf: f64,
}

// default clock options that correspond to simple linear regression
impl Default for ClockOptions {
  fn default() -> Self {
    Self {
      variance_factor: 0.0,
      variance_offset: 0.0,
      variance_offset_leaf: 1.0,
    }
  }
}

pub fn root_clock_model(graph: &ClockGraph) -> Result<ClockModel, Report> {
  // calculate a clock model from the root averages
  let root = graph.get_exactly_one_root()?;
  let root = root.read_arc().payload().read_arc();
  ClockModel::new(&root.total)
}

pub fn clock_regression_backward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_backward(|mut n| {
    let date = n.payload.date;
    // assemble the message that will be send to the parent
    let q_to_parent = if n.is_leaf {
      if n.payload.is_outlier {
        ClockSet::outlier_contribution()
      } else {
        ClockSet::leaf_contribution(date)
      }
    } else {
      let mut q_dest = ClockSet::default();
      for (c, e) in &n.children {
        q_dest += &e.read_arc().from_child;
      }
      q_dest
    };

    let is_leaf = n.is_leaf;
    let is_root = n.is_root;
    let edge_to_parent = n.get_exactly_one_parent_edge();

    if is_root {
      n.payload.total = q_to_parent;
    } else {
      // if not at the root, save the message to the parent on the edge
      let edge_to_parent = edge_to_parent.expect("Encountered a node without a parent edge");
      edge_to_parent.to_parent = q_to_parent.clone();
      let edge_len = edge_to_parent.weight().expect("Encountered an edge without a weight");
      let mut branch_variance = options.variance_factor * edge_len + options.variance_offset;
      // propagate the message to the parent along the edge (taking care of the speical case need for leafs)
      edge_to_parent.from_child = if is_leaf {
        branch_variance += options.variance_offset_leaf;
        ClockSet::leaf_contribution_to_parent(date, edge_len, branch_variance)
      } else {
        edge_to_parent.to_parent.propagate_averages(edge_len, branch_variance)
      };
    }
    return GraphTraversalContinuation::Continue;
  });
}

pub fn clock_regression_forward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_forward(|mut n| {
    if !n.is_root {
      // if not at the root, calculate the total message from the parent message, and the message sent to the parent
      let (parent, edge) = n.get_exactly_one_parent().unwrap();
      let edge = edge.read_arc();
      let edge_len = edge.weight().expect("Encountered an edge without a weight");
      let branch_variance = options.variance_factor * edge_len + options.variance_offset;

      let mut q_dest = edge.to_parent.clone();
      q_dest += edge.to_child.propagate_averages(edge_len, branch_variance);
      n.payload.total = q_dest.clone();
    }

    // place a message to the child onto each edge. this is the total message minus the message from the child
    for mut child_edge in n.child_edges {
      let mut q = n.payload.total.clone();
      q -= &child_edge.from_child;
      child_edge.to_child = q;
    }
    GraphTraversalContinuation::Continue
  });
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::clock::find_best_root::find_best_root;
  use crate::graph::node::Named;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::seq::div::{calculate_divs, OnlyLeaves};
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
