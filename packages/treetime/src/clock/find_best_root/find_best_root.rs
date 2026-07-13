use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::find_best_split::{FindRootResult, find_best_split};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::make_error;
use crate::payload::clock_set::ClockSet;
use crate::payload::traits::{ClockEdge, ClockNode};
use eyre::Report;
use log::{debug, info};
use rayon::prelude::*;
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::collections::container::get_exactly_one;

/// Find the best new root node
///
// Loop over all nodes, pick the one with the lowest chisq (and positive clock rate
// when force_positive is true), then optimize position along surrounding branches.
pub fn find_best_root<N, E, D>(
  graph: &Graph<N, E, D>,
  options: &ClockParams,
  params: &BranchPointOptimizationParams,
  force_positive: bool,
  objective: RootObjective,
) -> Result<FindRootResult, Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  info!("Starting root optimization with method: {params:?}, force_positive={force_positive}");

  let root = graph.get_exactly_one_root()?;
  let mut best_root_node = Arc::clone(&root);

  // Initialize with the current root, only accepting it if it has a positive clock rate
  // (or if force_positive is false)
  let root_clock_set = root.read_arc().payload().read_arc().clock_set().clone();
  let root_acceptable = !force_positive || has_positive_clock_rate(&root_clock_set);
  let mut best_chisq = if root_acceptable {
    objective.score(&root_clock_set)
  } else {
    f64::INFINITY
  };
  debug!(
    "Initial root chi-squared: {:.6e} (acceptable: {root_acceptable})",
    objective.score(&root_clock_set)
  );

  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    chisq: best_chisq,
    clock_set: root_clock_set,
  };

  // Find best node (with positive clock rate when force_positive is true)
  let mut node_count = 0;
  let mut improvements = 0;
  let mut rejected_negative_rate = 0;
  let candidates = graph
    .get_nodes()
    .par_iter()
    .map(|node| {
      let node_guard = node.read_arc();
      let payload_guard = node_guard.payload().read_arc();
      let clock_set = payload_guard.clock_set();
      let acceptable = !force_positive || has_positive_clock_rate(clock_set);
      (Arc::clone(node), acceptable.then(|| objective.score(clock_set)))
    })
    .collect::<Vec<_>>();
  for (n, score) in candidates {
    let Some(tmp_chisq) = score else {
      rejected_negative_rate += 1;
      node_count += 1;
      continue;
    };
    if tmp_chisq < best_chisq {
      improvements += 1;
      debug!("Found better node {improvements}: chi-squared improved from {best_chisq:.6e} to {tmp_chisq:.6e}");
      best_chisq = tmp_chisq;
      best_root_node = n;
    }
    node_count += 1;
  }
  debug!(
    "Evaluated {node_count} nodes, found {improvements} improvements, \
     rejected {rejected_negative_rate} with negative rate, best chi-squared: {best_chisq:.6e}"
  );
  let best_root_node = best_root_node.read_arc();

  // Check if some intermediate place on the parent branch is better
  if !best_root_node.is_root() {
    debug!("Optimizing position on parent branch");
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, *edge, options, params, objective)?;
    debug!(
      "Parent branch optimization result: chi-squared = {:.6e}, split = {:.6}",
      res.chisq, res.split
    );
    let split_acceptable = !force_positive || has_positive_clock_rate(&res.clock_set);
    if res.chisq < best_chisq && split_acceptable {
      debug!(
        "Parent branch optimization improved chi-squared from {:.6e} to {:.6e}",
        best_chisq, res.chisq
      );
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  // Check if some place on a child branch is better
  for (child_branch_count, e) in best_root_node.outbound().iter().enumerate() {
    debug!("Optimizing position on child branch {child_branch_count}");
    let res = find_best_split(graph, *e, options, params, objective)?;
    debug!(
      "Child branch {} optimization result: chi-squared = {:.6e}, split = {:.6}",
      child_branch_count, res.chisq, res.split
    );
    let split_acceptable = !force_positive || has_positive_clock_rate(&res.clock_set);
    if res.chisq < best_chisq && split_acceptable {
      debug!(
        "Child branch {} optimization improved chi-squared from {:.6e} to {:.6e}",
        child_branch_count, best_chisq, res.chisq
      );
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  if force_positive && !has_positive_clock_rate(&best_res.clock_set) {
    return make_error!(
      "Clock rate is negative for all root positions. \
       The data may lack temporal signal. Please specify --clock-rate explicitly."
    );
  }

  info!(
    "Root optimization completed. Final chi-squared: {:.6e}, split: {:.6}",
    best_res.chisq, best_res.split
  );

  Ok(best_res)
}

fn has_positive_clock_rate(clock_set: &ClockSet) -> bool {
  let det = clock_set.determinant();
  det > 0.0 && clock_set.clock_rate(det) > 0.0
}
