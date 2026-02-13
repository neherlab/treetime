use crate::commands::ancestral::marginal::update_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::{RerootChanges, RerootParams};
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
/// Performs clock-based rerooting, then calls `apply_reroot` on each partition
/// with bundled topology changes (edge split, edge merge, inverted edges).
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
  let clock_reroot_result =
    estimate_clock_model_with_reroot_policy(graph, clock_params, clock_rate, false, branch_params, &reroot_params)
      .wrap_err("Failed to estimate clock model with reroot")?;

  let new_root_key = graph.get_exactly_one_root()?.read_arc().key();

  // Update partitions if we have any and rerooting occurred
  if let Some(reroot_result) = &clock_reroot_result.reroot_result {
    if !partitions.is_empty() {
      // Build inverted edge keys: edges on path from old root to new root (if old root still exists)
      let inverted_edge_keys = if new_root_key != old_root_key && graph.get_node(old_root_key).is_some() {
        info!(
          "Root changed from {} to {} - computing reroot path",
          old_root_key.0, new_root_key.0
        );

        let path = graph
          .path_from_node_to_node(old_root_key, new_root_key)
          .wrap_err("Failed to compute path from old root to new root")?;

        path
          .iter()
          .filter_map(|(_, edge)| edge.as_ref().map(|e| e.read_arc().key()))
          .collect_vec()
      } else {
        vec![]
      };

      let changes = RerootChanges {
        edge_split: reroot_result.edge_split.clone(),
        edge_merge: reroot_result.edge_merge.clone(),
        inverted_edge_keys,
      };

      info!("Applying reroot changes to {} partitions", partitions.len());
      for partition in partitions {
        partition
          .write_arc()
          .apply_reroot(&changes)
          .wrap_err("Failed to apply reroot changes to partition")?;
      }

      // Recompute marginal messages after partition state update
      update_marginal(graph, partitions).wrap_err("Failed to update marginal after reroot")?;
    }
  }

  Ok(clock_reroot_result.clock_model)
}
