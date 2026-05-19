use crate::ancestral::marginal::update_marginal;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::timetree::EdgeTimetree;
use crate::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::reroot::RerootChanges;

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
  force_positive_rate: bool,
) -> Result<ClockModel, Report> {
  let reroot_params = RerootParams {
    force_positive_rate,
    ..RerootParams::default()
  };

  info!(
    "Reroot params: split_edge={}, remove_trivial_root={}, force_positive_rate={force_positive_rate}",
    reroot_params.split_edge, reroot_params.remove_trivial_root
  );

  // Perform clock-based rerooting
  let clock_reroot_result = estimate_clock_model_with_reroot_policy(
    graph,
    clock_params,
    clock_rate,
    false,
    branch_params,
    &reroot_params,
    None,
  )
  .wrap_err("Failed to estimate clock model with reroot")?;

  if let Some(reroot_result) = &clock_reroot_result.reroot_result {
    if !partitions.is_empty() {
      let changes = RerootChanges {
        edge_split: reroot_result.edge_split.clone(),
        edge_merge: reroot_result.edge_merge.clone(),
        inverted_edge_keys: reroot_result.inverted_edge_keys.clone(),
      };

      info!("Applying reroot changes to {} partitions", partitions.len());
      for partition in partitions {
        partition
          .write_arc()
          .apply_reroot(&changes)
          .wrap_err("Failed to apply reroot changes to partition")?;
      }

      update_marginal(graph, partitions).wrap_err("Failed to update marginal after reroot")?;
    }
  }

  Ok(clock_reroot_result.clock_model)
}
