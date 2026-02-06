use crate::commands::ancestral::marginal_unified::update_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::RerootParams;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::GraphTimetree;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;

/// Reroot tree for optimal temporal signal and update partition state.
///
/// Builds reroot policy from partition capabilities, performs clock-based rerooting,
/// then updates partition states along the reroot path.
pub fn reroot_tree(
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_params: &ClockParams,
  clock_rate: Option<f64>,
  branch_params: &BranchPointOptimizationParams,
) -> Result<ClockModel, Report> {
  let old_root_key = graph.get_exactly_one_root()?.read_arc().key();

  // Use default reroot params - always allow edge split and trivial root removal
  let reroot_params = RerootParams::default();

  info!(
    "Reroot params: split_edge={}, remove_trivial_root={}",
    reroot_params.split_edge, reroot_params.remove_trivial_root
  );

  // Perform clock-based rerooting
  let clock_model =
    estimate_clock_model_with_reroot_policy(graph, clock_params, clock_rate, false, branch_params, &reroot_params)
      .wrap_err("Failed to estimate clock model with reroot")?;

  let new_root_key = graph.get_exactly_one_root()?.read_arc().key();

  // If root changed and we have partitions, update partition state
  // Note: old_root may have been removed if it became trivial after reroot
  let old_root_exists = graph.get_node(old_root_key).is_some();
  if new_root_key != old_root_key && !partitions.is_empty() && old_root_exists {
    info!(
      "Root changed from {} to {} - updating partitions",
      old_root_key.0, new_root_key.0
    );

    // Compute path from old root to new root (post-reroot graph)
    // After reroot, new_root is at the top, so we walk FROM old_root TO new_root (upward)
    let path = graph
      .path_from_node_to_node(old_root_key, new_root_key)
      .wrap_err("Failed to compute path from old root to new root")?;

    // Convert to key-based representation for partition reroot
    let path_keys = path
      .iter()
      .map(|(node, edge)| {
        let node_key = node.read_arc().key();
        let edge_key = edge.as_ref().map(|e| e.read_arc().key());
        (node_key, edge_key)
      })
      .collect_vec();

    // Update each partition's state along the reroot path
    for partition in partitions {
      partition
        .write_arc()
        .reroot_partition_node_only(graph, old_root_key, new_root_key, &path_keys)
        .wrap_err("Failed to update partition after reroot")?;
    }

    // Recompute marginal messages after partition state update
    update_marginal(graph, partitions).wrap_err("Failed to update marginal after reroot")?;
  }

  Ok(clock_model)
}
