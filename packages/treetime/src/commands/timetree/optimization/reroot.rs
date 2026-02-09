use crate::commands::ancestral::marginal_unified::update_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::RerootParams;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::node::GraphNodeKey;
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
/// Performs clock-based rerooting, then updates partition states:
/// 1. Handle edge split if one occurred
/// 2. Handle edge merge if old root was removed
/// 3. Update partition state along reroot path
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
      // Handle edge split if one occurred
      if let Some(edge_split) = &reroot_result.edge_split {
        info!("Handling edge split: new node {}", edge_split.new_node_key.0);
        for partition in partitions {
          partition
            .write_arc()
            .handle_edge_split(edge_split)
            .wrap_err("Failed to handle edge split in partition")?;
        }
      }

      // Handle edge merge if old root was removed
      if let Some(edge_merge) = &reroot_result.edge_merge {
        info!("Handling edge merge: removed node {}", edge_merge.removed_node_key.0);
        for partition in partitions {
          partition
            .write_arc()
            .handle_edge_merge(edge_merge)
            .wrap_err("Failed to handle edge merge in partition")?;
        }
      }

      // Update partition state along reroot path if root changed and old root still exists
      let old_root_exists = graph.get_node(old_root_key).is_some();
      if new_root_key != old_root_key && old_root_exists {
        info!(
          "Root changed from {} to {} - updating partitions",
          old_root_key.0, new_root_key.0
        );

        // Compute path from old root to new root (post-reroot graph)
        let path = graph
          .path_from_node_to_node(old_root_key, new_root_key)
          .wrap_err("Failed to compute path from old root to new root")?;

        let path_keys = path
          .iter()
          .map(|(node, edge)| {
            let node_key = node.read_arc().key();
            let edge_key = edge.as_ref().map(|e| e.read_arc().key());
            (node_key, edge_key)
          })
          .collect_vec();

        for partition in partitions {
          partition
            .write_arc()
            .update_partition_after_reroot(graph, old_root_key, new_root_key, &path_keys)
            .wrap_err("Failed to update partition after reroot")?;
        }
      } else if new_root_key != old_root_key {
        // Old root was removed - use new root as reference point
        info!(
          "Root changed from {} (removed) to {} - updating partitions with new root only",
          old_root_key.0, new_root_key.0
        );

        // Path is just the new root node with no edges
        let path_keys: Vec<(GraphNodeKey, Option<GraphEdgeKey>)> = vec![(new_root_key, None)];

        for partition in partitions {
          partition
            .write_arc()
            .update_partition_after_reroot(graph, old_root_key, new_root_key, &path_keys)
            .wrap_err("Failed to update partition after reroot")?;
        }
      }

      // Recompute marginal messages after partition state update
      update_marginal(graph, partitions).wrap_err("Failed to update marginal after reroot")?;
    }
  }

  Ok(clock_reroot_result.clock_model)
}
