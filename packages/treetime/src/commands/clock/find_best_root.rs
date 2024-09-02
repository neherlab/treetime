use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{clock_regression_backward, clock_regression_forward, ClockOptions};
use crate::graph::edge::{GraphEdgeKey, Weighted};
use crate::graph::node::Named;
use crate::utils::container::get_exactly_one;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Debug, Serialize, Deserialize)]
pub struct FindRootResult {
  pub edge: Option<GraphEdgeKey>,

  /// The best root might be somewhere half-way on an existing edge.
  /// The position of this best root is the split. If this is 0 or 1, then we reroot on either the parent (source) or
  /// the child (target) of the edge. If this is 0 < x < 1, we put a new node at that point, and reroot on that new node.
  pub split: f64,

  pub clock: ClockModel,
}

/// Find the best new root node
pub fn find_best_root(graph: &ClockGraph, options: &ClockOptions) -> Result<FindRootResult, Report> {
  clock_regression_backward(graph, options);
  clock_regression_forward(graph, options);

  // Loop over all nodes, pick the one with the lowest chisq, then optimize position along surrounding branches.
  let root = graph.get_exactly_one_root()?;
  let mut best_root_node = Arc::clone(&root);

  let root = root.read_arc().payload().read_arc();
  let clock = ClockModel::new(&root.total)?;
  let mut best_chisq = clock.chisq();

  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    clock,
  };

  // find best node
  for n in graph.get_nodes() {
    let tmp_chisq = {
      let n = n.read_arc().payload().read_arc();
      ClockModel::new(&n.total)?.chisq() // FIXME: we might not need to calculate the whole model here, only the .chisq()
    };
    if tmp_chisq < best_chisq {
      best_chisq = tmp_chisq;
      best_root_node = Arc::clone(&n);
    }
  }

  let best_root_node = best_root_node.read_arc();

  // Check if someone on parent branch is better
  if !best_root_node.is_root() {
    // TODO: should we loop over parents as well? (just like over children in the loop just below)
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, *edge, options)?;
    if res.clock.chisq() < best_chisq {
      best_chisq = res.clock.chisq();
      best_res = res;
    }
  }

  // Check if someone on a child branch is better
  for e in best_root_node.outbound() {
    let res = find_best_split(graph, *e, options)?;
    if res.clock.chisq() < best_chisq {
      best_chisq = res.clock.chisq();
      best_res = res;
    }
  }

  Ok(best_res)
}

fn find_best_split(graph: &ClockGraph, edge: GraphEdgeKey, options: &ClockOptions) -> Result<FindRootResult, Report> {
  let edge = graph.get_edge(edge).expect("Edge not found");

  // Get clock data from both ends of the edge.
  let n = graph.get_node(edge.read_arc().target()).expect("Node not found");
  let n_clock = &n.read_arc().payload().read_arc().to_parent;
  let n = n.read_arc().payload().read_arc();
  let n_name = n.name().expect("Encountered node without a name").as_ref().to_owned();

  let p = graph.get_node(edge.read_arc().source()).expect("Node not found");
  let p_clock = &p.read_arc().payload().read_arc().to_children[&n_name];

  // Precalculate branch values for the edge.
  let branch_length = edge
    .read_arc()
    .payload()
    .read_arc()
    .weight()
    .expect("Encountered an edge without a weight");

  let branch_variance = options.variance_factor * branch_length + options.variance_offset;

  // Interrogate different positions along the branch
  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_clock = ClockModel::default();

  // TODO: arbitrary choice for now, should optimize
  for x in Array1::linspace(0.0, 1.0, 11) {
    let Q = p_clock.propagate_averages(branch_length * x, branch_variance * x)
      + n_clock.propagate_averages(branch_length * (1.0 - x), branch_variance * (1.0 - x));
    let clock_model = ClockModel::new(&Q)?;
    if clock_model.chisq() < best_chisq {
      best_chisq = clock_model.chisq();
      best_split = x;
      best_clock = clock_model;
    }
  }

  Ok(FindRootResult {
    edge: Some(edge.read_arc().key()),
    split: best_split,
    clock: best_clock,
  })
}
