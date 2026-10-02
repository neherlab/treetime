use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::clock_filter::{ClockFilterResult, clock_filter};
use crate::clock::clock_regression::{ClockFit, ClockVarianceParams};
use crate::clock::clock_state::ClockInputs;
use crate::clock::reroot::RerootParams;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts};
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::BranchOptMethod;
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::progress::{LogSink, StageSink};
use crate::progress_info;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::inference::time_inference::likely_times;
use crate::timetree::optimization::outliers::report_outliers;
use crate::timetree::optimization::reroot::{RerootedTree, reroot_tree};
use crate::timetree::pipeline::{TimetreeContext, TimetreeParams};
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

const TIMETREE_PRE_STEP_DAMPING: f64 = 0.75;

const PRE_LOOP_SCHEDULE: [PreLoopStep; 6] = [
  PreLoopStep::MlOptimizePreReroot,
  PreLoopStep::RerootPreAncestral,
  PreLoopStep::ClockFilter,
  PreLoopStep::MlOptimizePostReroot,
  PreLoopStep::InitialRoundCheckpoint,
  PreLoopStep::RerootPostAncestral,
];

pub(crate) fn run_pre_loop(
  inputs: &PreLoopInputs<'_>,
  state: PreLoopState,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<PreLoopState, Report> {
  PRE_LOOP_SCHEDULE.into_iter().try_fold(state, |state, step| {
    run_pre_loop_step(step, inputs, state, cancel, stages, log)
  })
}

#[derive(Clone, Copy)]
enum PreLoopStep {
  MlOptimizePreReroot,
  RerootPreAncestral,
  ClockFilter,
  MlOptimizePostReroot,
  InitialRoundCheckpoint,
  RerootPostAncestral,
}

pub(crate) struct PreLoopInputs<'a> {
  pub params: &'a TimetreeParams,
  pub context: &'a TimetreeContext,
  pub names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  pub has_alignment: bool,
}

pub(crate) struct PreLoopState {
  pub graph: Graph,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub branch_model: BranchModel,
  pub clock_fit: ClockFit,
  pub outliers: BTreeSet<GraphNodeKey>,
  pub filter_divergences: Option<BTreeMap<GraphNodeKey, f64>>,
}

impl PreLoopState {
  pub(crate) fn new(
    graph: Graph,
    branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
    branch_model: BranchModel,
    clock_fit: ClockFit,
  ) -> Self {
    Self {
      graph,
      branch_lengths,
      branch_model,
      clock_fit,
      outliers: BTreeSet::new(),
      filter_divergences: None,
    }
  }
}

fn run_pre_loop_step(
  step: PreLoopStep,
  inputs: &PreLoopInputs<'_>,
  state: PreLoopState,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<PreLoopState, Report> {
  let params = inputs.params;
  match step {
    PreLoopStep::MlOptimizePreReroot if inputs.has_alignment => ml_optimize(
      state,
      params.no_indels,
      "### ML branch-length optimization (pre-reroot)",
      "ML branch-length optimization (pre-reroot) failed",
      log,
    ),
    PreLoopStep::RerootPreAncestral if !params.keep_root => {
      progress_info!(log, "First reroot (pre-ancestral)");
      reroot(
        inputs,
        state,
        RerootOutliers::None,
        &ClockVarianceParams::default(),
        log,
      )
      .wrap_err("Failed to reroot tree (pre-ancestral)")
    },
    PreLoopStep::ClockFilter if params.clock_filter > 0.0 => filter_clock_outliers(inputs, state, log),
    PreLoopStep::MlOptimizePostReroot if inputs.has_alignment => {
      if matches!(state.branch_model, BranchModel::Input) {
        progress_info!(log, "Using input branch lengths for timetree inference");
        return Ok(state);
      }
      ml_optimize(
        state,
        params.no_indels,
        "### ML branch-length optimization (post-reroot)",
        "ML branch-length optimization (post-reroot) failed",
        log,
      )
    },
    PreLoopStep::InitialRoundCheckpoint => {
      cancel.check()?;
      stages.report("Initial timetree inference", 0.2, "");
      progress_info!(log, "### TreeTime: initial round");
      Ok(state)
    },
    PreLoopStep::RerootPostAncestral if !params.keep_root => {
      progress_info!(log, "Reroot (post-ancestral)");
      reroot(
        inputs,
        state,
        RerootOutliers::Filtered,
        &inputs.context.covariation_clock_params,
        log,
      )
      .wrap_err("Failed to reroot tree (post-ancestral)")
    },
    PreLoopStep::MlOptimizePreReroot
    | PreLoopStep::RerootPreAncestral
    | PreLoopStep::ClockFilter
    | PreLoopStep::MlOptimizePostReroot
    | PreLoopStep::RerootPostAncestral => Ok(state),
  }
}

fn ml_optimize(
  state: PreLoopState,
  no_indels: bool,
  banner: &str,
  failure: &'static str,
  log: &dyn LogSink,
) -> Result<PreLoopState, Report> {
  let BranchModel::Marginal(partition) = state.branch_model else {
    return Ok(state);
  };
  progress_info!(log, "{banner}");
  let (partition, _) = partition.marginal_update(&state.graph, &branch_lengths_or_zero(&state.branch_lengths))?;
  let (partition, branch_lengths) =
    optimize_branch_lengths(&state.graph, partition, state.branch_lengths, no_indels).wrap_err(failure)?;
  Ok(PreLoopState {
    branch_model: BranchModel::Marginal(partition),
    branch_lengths,
    ..state
  })
}

fn reroot(
  inputs: &PreLoopInputs<'_>,
  state: PreLoopState,
  outliers: RerootOutliers,
  clock_params: &ClockVarianceParams,
  log: &dyn LogSink,
) -> Result<PreLoopState, Report> {
  let params = inputs.params;
  let no_outliers = BTreeSet::new();
  let outliers = match outliers {
    RerootOutliers::None => &no_outliers,
    RerootOutliers::Filtered => &state.outliers,
  };
  let RerootedTree {
    graph,
    branch_lengths,
    branch_model,
    clock_fit,
  } = reroot_tree(
    state.graph,
    state.branch_lengths,
    state.branch_model,
    &inputs.context.date_constraints,
    outliers,
    clock_params,
    params.clock_rate,
    &inputs.context.branch_params,
    &RerootParams::new(params.reroot_spec.clone(), !params.allow_negative_rate),
    inputs.names,
    log,
  )?;
  Ok(PreLoopState {
    graph,
    branch_lengths,
    branch_model,
    clock_fit,
    ..state
  })
}

#[derive(Clone, Copy)]
enum RerootOutliers {
  None,
  Filtered,
}

fn filter_clock_outliers(
  inputs: &PreLoopInputs<'_>,
  state: PreLoopState,
  log: &dyn LogSink,
) -> Result<PreLoopState, Report> {
  let graph = &state.graph;
  let given_dates = likely_times(graph, &inputs.context.date_constraints, None)?;
  let clock_inputs = ClockInputs::from_times(graph, &given_dates, &BTreeMap::new());
  let ClockFilterResult {
    outliers,
    divergences,
    iqd,
  } = clock_filter(
    graph,
    &clock_inputs,
    &state.clock_fit.model,
    &state.branch_lengths,
    inputs.params.clock_filter,
    log,
  )?;
  report_outliers(
    graph,
    &outliers,
    &divergences,
    &state.clock_fit.model,
    iqd,
    &given_dates,
    inputs.names,
    log,
  );
  Ok(PreLoopState {
    outliers,
    filter_divergences: Some(divergences),
    ..state
  })
}

fn optimize_branch_lengths(
  graph: &Graph,
  partition: MarginalReconstruction,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  no_indels: bool,
) -> Result<(MarginalReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
  let old_branch_lengths = branch_lengths;
  let mut branch_lengths = old_branch_lengths.clone();

  let total_length = partition.sequence_length();
  let contributions = gather_edge_contributions(graph, &partition)?;
  let indel_counts = gather_edge_indel_counts(graph, &partition);
  if no_indels {
    run_optimize_mixed_inner(
      graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
      0.0,
      true,
      &mut branch_lengths,
    )
    .wrap_err("ML branch-length optimization pre-step failed")?;
  } else {
    run_optimize_mixed(
      graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
      &mut branch_lengths,
    )
    .wrap_err("ML branch-length optimization pre-step failed")?;
  }

  apply_damping(&mut branch_lengths, &old_branch_lengths, TIMETREE_PRE_STEP_DAMPING, 0);
  let (partition, _) = partition.marginal_update(graph, &branch_lengths_or_zero(&branch_lengths))?;
  Ok((partition, branch_lengths))
}
