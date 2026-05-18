use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::find_best_root::find_best_root;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::partition::algo::topology_cleanup::reroot::{
  self as topology_reroot, remove_node_if_trivial, split_edge,
};
use crate::partition::payload::clock_set::ClockSet;
use crate::partition::payload::traits::{ClockEdge, ClockNode};
use approx::ulps_eq;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

pub use topology_reroot::{EdgeMergeInfo, EdgeSplitInfo, RerootChanges, RerootResult};

/// Controls reroot behavior for graph topology changes.
#[derive(Clone, Debug, Serialize, Deserialize, SmartDefault)]
pub struct RerootParams {
  /// Allow creating a new node by splitting an edge during reroot.
  /// When false, reroot will snap to the nearest existing node endpoint.
  #[default = true]
  pub split_edge: bool,

  /// Remove the old root node if it becomes trivial (one parent, one child) after reroot.
  /// When false, the old root is preserved even if trivial.
  #[default = true]
  pub remove_trivial_root: bool,

  /// Only accept root positions with positive estimated clock rate.
  /// When false, the best chi-squared root is accepted regardless of rate sign.
  /// Use false for pre-filter steps where outliers may cause negative rates at all positions.
  #[default = true]
  pub force_positive_rate: bool,
}

pub fn reroot_in_place<N, E, D>(
  graph: &mut Graph<N, E, D>,
  options: &ClockParams,
  params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
) -> Result<RerootResult, Report>
where
  N: GraphNode + ClockNode + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  let FindRootResult {
    edge, split, clock_set, ..
  } = find_best_root(graph, options, params, reroot_params.force_positive_rate)?;

  let old_root_key = { graph.get_exactly_one_root()?.read_arc().key() };
  let Some(edge_key) = edge else {
    // Already at the best root
    return Ok(RerootResult {
      new_root_key: old_root_key,
      edge_split: None,
      edge_merge: None,
    });
  };

  // Extract edge endpoints before the edge is removed by split_edge
  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    (edge.read_arc().source(), edge.read_arc().target())
  };

  // Determine where to place the new root based on the split position along the edge
  // - split = 0.0 means root at target
  // - split = 1.0 means root at source
  // - 0.0 < split < 1.0 means somewhere in the middle - i.e. create a new node
  let (new_root_key, edge_split) = if ulps_eq!(split, 0.0, max_ulps = 5) {
    (target_key, None)
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    (source_key, None)
  } else if reroot_params.split_edge {
    let split_info = create_new_root_node(graph, edge_key, split, clock_set)?;
    (split_info.new_node_key, Some(split_info))
  } else {
    // Edge split disallowed: snap to nearest endpoint
    (if split < 0.5 { target_key } else { source_key }, None)
  };

  let edge_merge = if new_root_key != old_root_key {
    apply_reroot(graph, old_root_key, new_root_key, options)?;

    if reroot_params.remove_trivial_root {
      remove_node_if_trivial(graph, old_root_key)?
    } else {
      None
    }
  } else {
    None
  };

  Ok(RerootResult {
    new_root_key,
    edge_split,
    edge_merge,
  })
}

/// Create new root node by splitting the edge into two, then setting clock data on the new node.
fn create_new_root_node<N, E, D>(
  graph: &mut Graph<N, E, D>,
  edge_key: GraphEdgeKey,
  split: f64,
  clock_set: ClockSet,
) -> Result<EdgeSplitInfo, Report>
where
  N: GraphNode + ClockNode + Default,
  E: GraphEdge + ClockEdge + Default,
  D: Send + Sync,
{
  let split_info = split_edge(graph, edge_key, split)?;

  // Set clock-specific data on the new split node
  let new_node = graph
    .get_node(split_info.new_node_key)
    .expect("New split node not found");
  new_node.write_arc().payload().write_arc().set_clock_set(clock_set);

  Ok(split_info)
}

/// Modify graph topology to make the newly identified root the actual root,
/// then update clock-specific edge messages.
fn apply_reroot<N, E, D>(
  graph: &mut Graph<N, E, D>,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
  options: &ClockParams,
) -> Result<(), Report>
where
  N: GraphNode + ClockNode,
  E: GraphEdge + ClockEdge,
  D: Send + Sync,
{
  // Invert edges along the path (generic topology operation)
  let inverted_edge_keys = topology_reroot::apply_reroot_topology(graph, old_root_key, new_root_key)?;

  // Update clock-specific edge messages for inverted edges
  for edge_key in &inverted_edge_keys {
    let edge = graph.get_edge(*edge_key).expect("Inverted edge not found");
    let mut edge_payload = edge.read_arc().payload().write_arc();
    let edge_len = edge_payload.branch_length().unwrap();
    let branch_variance = options.variance_factor * edge_len + options.variance_offset;
    let tmp_to_parent = edge_payload.to_parent().clone();
    *edge_payload.to_parent_mut() = edge_payload.to_child().clone();
    *edge_payload.to_child_mut() = tmp_to_parent;
    *edge_payload.from_child_mut() = edge_payload.to_parent().propagate_averages(edge_len, branch_variance);
  }

  Ok(())
}
