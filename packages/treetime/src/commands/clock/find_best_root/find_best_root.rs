use crate::commands::clock::clock_regression::ClockParams;
use crate::commands::clock::clock_traits::{ClockEdge, ClockNode};
use crate::commands::clock::find_best_root::find_best_split::{FindRootResult, find_best_split};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use eyre::Report;
use log::{debug, info};
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::container::get_exactly_one;

/// Find the best new root node
///
// Loop over all nodes, pick the one with the lowest chisq, then optimize position along surrounding branches.
pub fn find_best_root<N, E, D>(
  graph: &Graph<N, E, D>,
  options: &ClockParams,
  params: &BranchPointOptimizationParams,
) -> Result<FindRootResult, Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  info!("Starting root optimization with method: {params:?}");

  let root = graph.get_exactly_one_root()?;
  let mut best_root_node = Arc::clone(&root);

  // Initialize with the current root
  let root = root.read_arc().payload().read_arc();
  let mut best_chisq = root.clock_set().chisq();
  debug!("Initial root chi-squared: {best_chisq:.6e}");

  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    chisq: best_chisq,
    clock_set: root.clock_set().clone(),
  };

  // Find best node
  let mut node_count = 0;
  let mut improvements = 0;
  for n in graph.get_nodes() {
    let tmp_chisq = n.read_arc().payload().read_arc().clock_set().chisq();
    if tmp_chisq < best_chisq {
      improvements += 1;
      debug!("Found better node {improvements}: chi-squared improved from {best_chisq:.6e} to {tmp_chisq:.6e}");
      best_chisq = tmp_chisq;
      best_root_node = Arc::clone(&n);
    }
    node_count += 1;
  }
  debug!("Evaluated {node_count} nodes, found {improvements} improvements, best chi-squared: {best_chisq:.6e}");
  let best_root_node = best_root_node.read_arc();

  // Check if some intermediate place on the parent branch is better
  if !best_root_node.is_root() {
    debug!("Optimizing position on parent branch");
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, *edge, options, params)?;
    debug!(
      "Parent branch optimization result: chi-squared = {:.6e}, split = {:.6}",
      res.chisq, res.split
    );
    if res.chisq < best_chisq {
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
    let res = find_best_split(graph, *e, options, params)?;
    debug!(
      "Child branch {} optimization result: chi-squared = {:.6e}, split = {:.6}",
      child_branch_count, res.chisq, res.split
    );
    if res.chisq < best_chisq {
      debug!(
        "Child branch {} optimization improved chi-squared from {:.6e} to {:.6e}",
        child_branch_count, best_chisq, res.chisq
      );
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  info!(
    "Root optimization completed. Final chi-squared: {:.6e}, split: {:.6}",
    best_res.chisq, best_res.split
  );

  Ok(best_res)
}
