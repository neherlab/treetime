use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_regression::{clock_regression_backward, clock_regression_forward, ClockOptions};
use crate::graph::edge::{GraphEdgeKey, Weighted};
use crate::utils::container::get_exactly_one;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

use super::clock_set::ClockSet;

#[derive(Debug, Serialize, Deserialize)]
pub struct FindRootResult {
  pub edge: Option<GraphEdgeKey>,

  /// The best root might be somewhere half-way on an existing edge.
  /// The position of this best root is the split. If this is 0 or 1, then we reroot on either the parent (source) or
  /// the child (target) of the edge. If this is 0 < x < 1, we put a new node at that point, and reroot on that new node.
  pub split: f64,

  pub total: ClockSet,

  pub chisq: f64,
}

/// Find the best new root node
pub fn find_best_root(graph: &ClockGraph, options: &ClockOptions) -> Result<FindRootResult, Report> {
  // Loop over all nodes, pick the one with the lowest chisq, then optimize position along surrounding branches.
  let root = graph.get_exactly_one_root()?;
  let mut best_root_node = Arc::clone(&root);

  // inialize with the current root
  let root = root.read_arc().payload().read_arc();
  let mut best_chisq = root.total.chisq();
  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    chisq: best_chisq,
    total: root.total.clone(),
  };

  // find best node
  for n in graph.get_nodes() {
    let tmp_chisq = n.read_arc().payload().read_arc().total.chisq();
    if tmp_chisq < best_chisq {
      best_chisq = tmp_chisq;
      best_root_node = Arc::clone(&n);
    }
  }
  let best_root_node = best_root_node.read_arc();

  // Check if some intermediate place on the parent branch is better
  if !best_root_node.is_root() {
    // TODO: should we loop over parents as well? (just like over children in the loop just below)
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, *edge, options)?;
    if res.chisq < best_chisq {
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  // Check if some place on a child branch is better
  for e in best_root_node.outbound() {
    let res = find_best_split(graph, *e, options)?;
    if res.chisq < best_chisq {
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  Ok(best_res)
}

fn find_best_split(graph: &ClockGraph, edge: GraphEdgeKey, options: &ClockOptions) -> Result<FindRootResult, Report> {
  let edge = graph.get_edge(edge).expect("Edge not found");
  let edge_payload = edge.read_arc().payload().read_arc();

  // check whether the target node is a leaf
  let is_terminal = graph.get_node(edge.read_arc().target()).unwrap().read_arc().is_leaf();

  // Precalculate branch values for the edge.
  let branch_length = edge_payload.weight().expect("Encountered an edge without a weight");

  let branch_variance = options.variance_factor * branch_length + options.variance_offset;

  // Interrogate different positions along the branch
  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_total: ClockSet = ClockSet::default();

  // TODO: arbitrary choice for now, should optimize
  for x in Array1::linspace(0.0, 1.0, 11) {
    // determine contribution of child/target first -- terminal nodes need special handling
    let child_contribution = if is_terminal {
      // hack: we could get this off the node as well
      let date = Some(edge_payload.to_parent.t_sum());
      ClockSet::leaf_contribution_to_parent(
        date,
        branch_length * (1.0 - x),
        branch_variance * (1.0 - x) + options.variance_offset_leaf,
      )
    } else {
      edge_payload
        .to_parent
        .propagate_averages(branch_length * (1.0 - x), branch_variance * (1.0 - x))
    };
    let clock_total = edge_payload
      .to_child
      .propagate_averages(branch_length * x, branch_variance * x)
      + child_contribution;
    let chisq = clock_total.chisq();
    if chisq < best_chisq {
      best_split = x;
      best_total = clock_total;
      best_chisq = chisq;
    }
  }

  Ok(FindRootResult {
    edge: Some(edge.read_arc().key()),
    split: best_split,
    total: best_total,
    chisq: best_chisq,
  })
}
