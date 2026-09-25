use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::clock::clock_regression::{ClockFit, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::partition::timetree::marginal::marginal_update_timetree;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::progress::ProgressSink;
use crate::progress_info;
use crate::timetree::inference::time_inference::likely_times;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
pub(crate) fn reroot_tree(
  graph: &mut Graph,
  constraints: &DateConstraints,
  clock_state: &mut ClockState,
  mut partitions: Vec<PartitionTimetree>,
  clock_params: &ClockVarianceParams,
  clock_rate: Option<f64>,
  branch_params: &BranchPointOptimizationParams,
  reroot_spec: &RerootSpec,
  force_positive_rate: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
) -> Result<(ClockFit, Vec<PartitionTimetree>), Report> {
  let reroot_params = RerootParams {
    spec: reroot_spec.clone(),
    force_positive_rate,
    ..RerootParams::default()
  };

  progress_info!(
    progress,
    "Reroot params: split_edge={}, remove_trivial_root={}, force_positive_rate={force_positive_rate}",
    reroot_params.split_edge,
    reroot_params.remove_trivial_root
  );

  clock_state.reseed_transitional(graph);
  let mut clock_inputs = ClockInputs::seed_from_times(graph, &likely_times(graph, constraints, None)?);
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
    progress,
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

      progress_info!(progress, "Applying reroot changes to {} partitions", partitions.len());
      partitions = partitions
        .into_iter()
        .map(|partition| partition.apply_reroot(&changes))
        .collect::<Result<Vec<_>, Report>>()
        .wrap_err("Failed to apply reroot changes to partition")?;

      (partitions, _) = marginal_update_timetree(graph, &branch_lengths_or_zero(branch_lengths), partitions)
        .wrap_err("Failed to update marginal after reroot")?;
    }
  }

  Ok((clock_reroot_result.into_clock_fit()?, partitions))
}
