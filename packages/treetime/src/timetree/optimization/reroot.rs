use crate::ancestral::marginal::profile_branch_lengths;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::partition::timetree::marginal::marginal_update_timetree;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::timetree_state::TimetreeState;
use eyre::{Report, WrapErr};
use log::info;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;

/// Reroot tree for optimal temporal signal and update partition state.
///
/// Performs clock-based rerooting, then calls `apply_reroot` on each partition
/// with bundled topology changes (edge split, edge merge, inverted edges).
#[allow(clippy::too_many_arguments)]
pub fn reroot_tree(
  graph: &mut Graph,
  constraints: &DateConstraints,
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
  // state into fresh clock inputs while preserving the value-resident divergence and outlier flag on
  // the clock results, so the regression excludes the leaves the clock filter marked. The reroot
  // rebuilds the returned results and remaps the inputs to match the new topology.
  clock_state.reseed_transitional(graph);
  let mut clock_inputs = ClockInputs::seed_from_times(graph, &timetree_state.likely_times(constraints));
  let (new_clock_state, clock_reroot_result) = estimate_clock_model_with_reroot_policy(
    graph,
    &mut clock_inputs,
    std::mem::take(clock_state),
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
  *clock_state = new_clock_state;

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

      marginal_update_timetree(graph, &profile_branch_lengths(branch_lengths), partitions)
        .wrap_err("Failed to update marginal after reroot")?;
    }
  }

  clock_reroot_result.into_clock_model()
}
