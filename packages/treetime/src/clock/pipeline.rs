use crate::cancel::Cancel;
use crate::clock::assign_dates::assign_dates;
use crate::clock::clock_filter::clock_filter_inplace;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockRegressionResult, gather_clock_regression_results};
use crate::progress::ProgressSink;
use eyre::{Report, WrapErr};
use log::info;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::date::DatesMap;

pub struct ClockPipelineParams {
  pub clock_params: ClockParams,
  pub clock_filter: f64,
  pub keep_root: bool,
  pub allow_negative_rate: bool,
  pub branch_params: BranchPointOptimizationParams,
  pub reroot_spec: RerootSpec,
}

pub struct ClockInput {
  pub graph: Graph,
  pub dates: DatesMap,
  /// Raw per-edge branch lengths captured from the Newick parse, keyed by edge id. Routed through
  /// estimation (the reroot search updates it in place) and read by the regression gather.
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct ClockOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub inputs: ClockInputs,
  #[serde(skip)]
  pub state: ClockState,
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

pub fn run(
  params: &ClockPipelineParams,
  mut input: ClockInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<ClockOutput, Report> {
  cancel.check()?;
  progress.report("Assigning dates", 0.1, "");
  let mut inputs = ClockInputs::new(&input.graph);
  assign_dates(&input.graph, &input.dates, &mut inputs, names)?;
  let state = ClockState::new(&input.graph);

  cancel.check()?;
  progress.report("Clock regression", 0.3, "");
  let mut branch_lengths = input.branch_lengths;
  let (mut state, clock_model, new_outliers) = estimate_clock_model_with_prefilter(
    &mut input.graph,
    &mut inputs,
    state,
    &params.clock_params,
    params.keep_root,
    &params.branch_params,
    params.clock_filter,
    params.allow_negative_rate,
    &params.reroot_spec,
    &mut branch_lengths,
    names,
  )?;

  if let Some(delta) = new_outliers {
    info!("Clock filter changed outlier status for {delta} leaf nodes");
  }

  // Estimation rerooted the tree in place: the best-root search updates `branch_lengths` for the new
  // edge set (a split adds two edges and drops one; merging the old trivial root drops two and adds
  // one), so the routed map is already current. Reproject the routed names onto the post-reroot node
  // set: a surviving node keeps its name (no clock pass renames), and the `N::default` node the split
  // adds is absent from the map and resolves to `None`. Clock, unlike timetree, does not name the new
  // root, so the map must not synthesize a name for it.
  let names: BTreeMap<GraphNodeKey, Option<String>> = input
    .graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, names.get(&key).cloned().flatten())
    })
    .collect();
  let regression_results =
    gather_clock_regression_results(&input.graph, &inputs, &mut state, &clock_model, &names, &branch_lengths)?;

  progress.report("Done", 1.0, "");
  Ok(ClockOutput {
    graph: input.graph,
    inputs,
    state,
    clock_model,
    regression_results,
    names,
    branch_lengths,
  })
}

#[allow(clippy::too_many_arguments)]
fn estimate_clock_model_with_prefilter(
  graph: &mut Graph,
  inputs: &mut ClockInputs,
  mut state: ClockState,
  options: &ClockParams,
  keep_root: bool,
  branch_params: &BranchPointOptimizationParams,
  clock_filter_threshold: f64,
  allow_negative_rate: bool,
  reroot_spec: &RerootSpec,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(ClockState, ClockModel, Option<i32>), Report> {
  let mut delta = None;
  if clock_filter_threshold > 0.0 {
    // Allow negative rates during pre-filter root finding. Some datasets (e.g. dengue/100)
    // have negative estimated rate at ALL root positions when outliers are included.
    // The pre-filter clock model only needs to be good enough for IQD-based outlier detection.
    let reroot_params = RerootParams {
      spec: reroot_spec.clone(),
      force_positive_rate: false,
      ..RerootParams::default()
    };
    let (new_state, result) = estimate_clock_model_with_reroot_policy(
      graph,
      inputs,
      state,
      options,
      None,
      keep_root,
      branch_params,
      &reroot_params,
      branch_lengths,
      None,
      names,
    )?;
    state = new_state;
    let regression = result.regression();
    if regression.clock_rate() < 0.0 {
      // IQD-based filtering uses |deviation| > IQD * threshold, so the absolute-value
      // comparison is slope-sign-invariant: outliers are identified by distance from the
      // fitted line regardless of slope direction.
      log::warn!(
        "Pre-filter clock rate is negative ({:.6e}). Outlier detection proceeds with this model.",
        regression.clock_rate()
      );
    }
    delta = Some(
      clock_filter_inplace(
        graph,
        inputs,
        &mut state,
        regression,
        branch_lengths,
        clock_filter_threshold,
      )?
      .new_outliers,
    );
  }

  let reroot_params = RerootParams {
    spec: reroot_spec.clone(),
    force_positive_rate: !allow_negative_rate,
    ..RerootParams::default()
  };
  let (state, result) = estimate_clock_model_with_reroot_policy(graph, inputs, state, options, None, keep_root, branch_params, &reroot_params, branch_lengths, None, names)
    .wrap_err_with(|| {
      if delta.is_some() {
        "Clock model estimation failed after outlier filtering. The pre-filter step removed outliers but the clock rate remains negative at all root positions.".to_owned()
      } else {
        "Clock model estimation failed".to_owned()
      }
    })?;
  // The clock command reports the root-to-tip regression without performing time
  // inference, so a non-positive rate is a valid (if temporally uninformative) result:
  // warn and continue rather than error. This matches v0 and makes `--allow-negative-rate`
  // and `--keep-root` usable. Timetree, which needs `time = div / rate`, still errors on a
  // non-positive rate (kb/decisions/timetree-rejects-negative-clock-rate.md).
  Ok((state, result.into_clock_model_allow_negative(), delta))
}
