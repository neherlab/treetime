use crate::ancestral::marginal::{marginal_update, profile_branch_lengths};
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::ClockState;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetree};
use crate::partition::traits::PartitionRerootOps;
use crate::timetree::timetree_state::TimetreeState;
use eyre::{Report, WrapErr};
use log::info;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;

/// Reroot tree for optimal temporal signal and update partition state.
///
/// Performs clock-based rerooting, then calls `apply_reroot` on each partition
/// with bundled topology changes (edge split, edge merge, inverted edges).
pub fn reroot_tree(
  graph: &mut GraphTimetree,
  clock_state: &mut ClockState,
  timetree_state: &TimetreeState,
  partitions: &mut [PartitionTimetree],
  clock_params: &ClockParams,
  clock_rate: Option<f64>,
  branch_params: &BranchPointOptimizationParams,
  reroot_spec: &RerootSpec,
  force_positive_rate: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<ClockModel, Report> {
  let reroot_params = RerootParams {
    spec: reroot_spec.clone(),
    force_positive_rate,
    ..RerootParams::default()
  };

  info!(
    "Reroot params: split_edge={}, remove_trivial_root={}, force_positive_rate={force_positive_rate}",
    reroot_params.split_edge, reroot_params.remove_trivial_root
  );

  // Perform clock-based rerooting on the threaded clock state. Re-read the node dates from the date
  // state while preserving the value-resident divergence and outlier flag, so the regression excludes
  // the leaves the clock filter marked. The reroot mutates the state's node and edge maps in place to
  // match the new topology.
  clock_state.reseed_transitional_from_times(graph, &timetree_state.likely_times(), &BTreeMap::new());
  let clock_reroot_result = estimate_clock_model_with_reroot_policy(
    graph,
    clock_state,
    clock_params,
    clock_rate,
    false,
    branch_params,
    &reroot_params,
    branch_lengths,
    None,
    names,
  )
  .wrap_err("Failed to estimate clock model with reroot")?;

  if let Some(reroot_result) = clock_reroot_result.reroot_result() {
    if !partitions.is_empty() {
      let changes = RerootChanges {
        edge_split: reroot_result.edge_split.clone(),
        edge_merge: reroot_result.edge_merge.clone(),
        inverted_edge_keys: reroot_result.inverted_edge_keys.clone(),
      };

      info!("Applying reroot changes to {} partitions", partitions.len());
      for partition in partitions.iter_mut() {
        partition
          .apply_reroot(&changes)
          .wrap_err("Failed to apply reroot changes to partition")?;
      }

      marginal_update(graph, &profile_branch_lengths(branch_lengths), partitions)
        .wrap_err("Failed to update marginal after reroot")?;
    }
  }

  clock_reroot_result.into_clock_model()
}
