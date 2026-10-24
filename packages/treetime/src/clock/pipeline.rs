use crate::branch_lengths::branch_length_or_zero;
use crate::cancel::Cancel;
use crate::clock::assign_dates::assign_dates;
use crate::clock::clock_filter::clock_filter;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockFitParams, ClockFitResult, ClockTree, fit_clock_model};
use crate::clock::clock_state::ClockInputs;
use crate::clock::divergence::root_to_node_divergences;
use crate::clock::rtt::{ClockRegressionResult, gather_clock_regression_results};
use crate::error::{NoPositiveRateRoot, OperationError};
use crate::progress::{LogSink, StageSink};
use crate::progress_warn;
use crate::reroot::date_free::{LengthsAdvice, reroot_date_free};
use crate::reroot::driver::{reroot_at_best_root, update_clock_inputs};
use crate::reroot::params::{RerootParams, RerootSpec};
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::date::DateConstraint;
use treetime_utils::make_internal_error;

pub fn run(
  params: &ClockParams,
  input: ClockInput,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<ClockOutput, OperationError> {
  cancel.check().map_err(OperationError::classify)?;
  stages.report("Assigning dates", 0.1, "");
  let mut inputs = ClockInputs::new(&input.graph);
  assign_dates(&input.graph, &input.dates, &mut inputs).map_err(OperationError::InvalidInput)?;

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Clock regression", 0.3, "");
  let mut names = input.names;
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
  ) = estimate_clock_model_with_prefilter(tree, params, &mut names, log).map_err(OperationError::classify)?;

  let outliers = filter_outliers.unwrap_or_default();

  let divergences = root_to_node_divergences(&graph, |edge_key| branch_length_or_zero(&branch_lengths, edge_key))
    .map_err(OperationError::classify)?;
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

pub struct ClockInput {
  pub graph: Graph,
  pub names: NodeNames,
  pub dates: BTreeMap<GraphNodeKey, DateConstraint>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug)]
pub struct ClockOutput {
  pub graph: Graph,
  pub inputs: ClockInputs,
  pub divergences: BTreeMap<GraphNodeKey, f64>,
  pub outliers: BTreeSet<GraphNodeKey>,
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
  pub names: NodeNames,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

fn estimate_clock_model_with_prefilter(
  tree: ClockTree,
  params: &ClockParams,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockModel, Option<BTreeSet<GraphNodeKey>>), Report> {
  let (clock_filter_threshold, allow_negative_rate) = (params.clock_filter, params.allow_negative_rate);
  let mut tree = tree;
  if let Some(spec) = params.reroot_spec.as_ref().filter(|spec| !spec.uses_dates()) {
    reroot_clock_tree_date_free(&mut tree, spec, names, log)?;
  }
  let reroot_spec = params.reroot_spec.as_ref().filter(|spec| spec.uses_dates());
  let filter_outliers = if clock_filter_threshold > 0.0 {
    let result = reroot_and_fit(&mut tree, reroot_spec, &BTreeSet::new(), true, names, log).map_err(|report| {
      if report.downcast_ref::<NoPositiveRateRoot>().is_some() {
        report.wrap_err(
          "The clock fit of the outlier filter requires a positive clock rate, also with --allow-negative-rate. Run \
           with --clock-filter=0 to skip the outlier filter",
        )
      } else {
        report
      }
    })?;
    let Some(regression) = result.regression() else {
      return make_internal_error!(
        "The clock fit before outlier filtering estimates the rate, so it must hold a regression"
      );
    };
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
    Some(filtered.outliers)
  } else {
    None
  };

  let no_outliers = BTreeSet::new();
  let outliers = filter_outliers.as_ref().unwrap_or(&no_outliers);
  let failure = if filter_outliers.is_some() {
    "Clock model estimation failed after outlier filtering. The pre-filter step removed outliers but the clock rate remains negative at all root positions."
  } else {
    "Clock model estimation failed"
  };
  let result = reroot_and_fit(&mut tree, reroot_spec, outliers, !allow_negative_rate, names, log).wrap_err(failure)?;
  Ok((tree, result.into_clock_model_allow_negative(log), filter_outliers))
}

fn reroot_clock_tree_date_free(
  tree: &mut ClockTree,
  spec: &RerootSpec,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let ClockTree {
    graph,
    branch_lengths,
    inputs,
  } = tree;
  let result = reroot_date_free(
    graph,
    branch_lengths,
    spec,
    &|key| inputs.likely_time(key).is_some(),
    LengthsAdvice::OptimizeTheTree,
    names,
    log,
  )?;
  update_clock_inputs(inputs, &result);
  Ok(())
}

fn reroot_and_fit(
  tree: &mut ClockTree,
  reroot_spec: Option<&RerootSpec>,
  outliers: &BTreeSet<GraphNodeKey>,
  force_positive_rate: bool,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<ClockFitResult, Report> {
  if let Some(spec) = reroot_spec {
    let ClockTree {
      graph,
      branch_lengths,
      inputs,
    } = tree;
    reroot_at_best_root(
      graph,
      branch_lengths,
      inputs,
      outliers,
      &RerootParams::new(spec.clone(), force_positive_rate),
      names,
      log,
    )?;
  }
  fit_clock_model(
    tree,
    &ClockFitParams {
      outliers,
      clock_rate: None,
      prev_clock_rate: None,
    },
    log,
  )
}

pub struct ClockParams {
  pub clock_filter: f64,
  pub allow_negative_rate: bool,
  pub reroot_spec: Option<RerootSpec>,
}
