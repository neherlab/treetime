use crate::branch_lengths::branch_length_or_zero;
use crate::cancel::Cancel;
use crate::clock::assign_dates::assign_dates;
use crate::clock::clock_filter::clock_filter;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::ClockInputs;
use crate::clock::divergence::root_to_node_divergences;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockRegressionResult, gather_clock_regression_results};
use crate::error::OperationError;
use crate::progress::{LogSink, StageSink};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::restrict_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::date::DatesMap;

pub fn run(
  params: &ClockParams,
  input: ClockInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<ClockOutput, OperationError> {
  cancel.check().map_err(OperationError::from_inference)?;
  stages.report("Assigning dates", 0.1, "");
  let mut inputs = ClockInputs::new(&input.graph);
  assign_dates(&input.graph, &input.dates, &mut inputs, names).map_err(OperationError::InvalidInput)?;

  cancel.check().map_err(OperationError::from_inference)?;
  stages.report("Clock regression", 0.3, "");
  let tree = ClockTree {
    graph: input.graph,
    branch_lengths: input.branch_lengths,
    inputs,
  };
  let (
    ClockTree {
      graph,
      branch_lengths,
      inputs,
    },
    clock_model,
    filter_outliers,
  ) = estimate_clock_model_with_prefilter(
    tree,
    &params.clock_params,
    params.keep_root,
    &params.branch_params,
    params.clock_filter,
    params.allow_negative_rate,
    &params.reroot_spec,
    names,
    log,
  )
  .map_err(OperationError::InvalidInput)?;

  if let Some(outliers) = &filter_outliers {
    progress_info!(log, "Clock filter flagged {} leaf nodes as outliers", outliers.len());
  }
  let outliers = filter_outliers.unwrap_or_default();

  let names = restrict_node_names(names, &graph);
  let divergences = root_to_node_divergences(&graph, |edge_key| branch_length_or_zero(&branch_lengths, edge_key))
    .map_err(OperationError::from_inference)?;
  let regression_results =
    gather_clock_regression_results(&graph, &inputs, &divergences, &outliers, &clock_model, &names);

  Ok(ClockOutput {
    graph,
    inputs,
    divergences,
    outliers,
    clock_model,
    regression_results,
    names,
    branch_lengths,
  })
}

pub struct ClockParams {
  pub clock_params: ClockVarianceParams,
  pub clock_filter: f64,
  pub keep_root: bool,
  pub allow_negative_rate: bool,
  pub branch_params: BranchPointOptimizationParams,
  pub reroot_spec: RerootSpec,
}

pub struct ClockInput {
  pub graph: Graph,
  pub dates: DatesMap,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct ClockOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub inputs: ClockInputs,
  #[serde(skip)]
  pub divergences: BTreeMap<GraphNodeKey, f64>,
  #[serde(skip)]
  pub outliers: BTreeSet<GraphNodeKey>,
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
fn estimate_clock_model_with_prefilter(
  tree: ClockTree,
  options: &ClockVarianceParams,
  keep_root: bool,
  branch_params: &BranchPointOptimizationParams,
  clock_filter_threshold: f64,
  allow_negative_rate: bool,
  reroot_spec: &RerootSpec,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockModel, Option<BTreeSet<GraphNodeKey>>), Report> {
  let (tree, filter_outliers) = if clock_filter_threshold > 0.0 {
    let reroot_params = RerootParams::new(reroot_spec.clone(), false);
    let (tree, result) = estimate_clock_model_with_reroot_policy(
      tree,
      &BTreeSet::new(),
      options,
      None,
      keep_root,
      branch_params,
      &reroot_params,
      None,
      names,
      log,
    )?;
    let regression = result.regression();
    if regression.clock_rate() < 0.0 {
      progress_warn!(
        log,
        "Pre-filter clock rate is negative ({:.6e}). Outlier detection proceeds with this model.",
        regression.clock_rate()
      );
    }
    let filtered = clock_filter(
      &tree.graph,
      &tree.inputs,
      regression,
      &tree.branch_lengths,
      clock_filter_threshold,
      log,
    )?;
    (tree, Some(filtered.outliers))
  } else {
    (tree, None)
  };

  let reroot_params = RerootParams::new(reroot_spec.clone(), !allow_negative_rate);
  let no_outliers = BTreeSet::new();
  let outliers = filter_outliers.as_ref().unwrap_or(&no_outliers);
  let failure = if filter_outliers.is_some() {
    "Clock model estimation failed after outlier filtering. The pre-filter step removed outliers but the clock rate remains negative at all root positions."
  } else {
    "Clock model estimation failed"
  };
  let (tree, result) = estimate_clock_model_with_reroot_policy(
    tree,
    outliers,
    options,
    None,
    keep_root,
    branch_params,
    &reroot_params,
    None,
    names,
    log,
  )
  .wrap_err(failure)?;
  Ok((tree, result.into_clock_model_allow_negative(log), filter_outliers))
}
