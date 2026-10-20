use crate::clock::clock_regression::ClockFit;
use crate::clock::gaussian::{GaussianClock, GaussianTree};
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::coalescent::skyline::SkylineParams;
use crate::progress::LogSink;
use crate::timetree::coalescent_timescale::{
  CoalescentMode, CoalescentSetup, LineageModel, coalescent_mode, coalescent_timescale,
};
use crate::timetree::inference::bad_branches::bad_leaves;
use crate::timetree::inference::result::{TimeInference, given_times};
use crate::timetree::inference::runner::{CLOCK_BRANCH_LENGTH_UNDAMPED, run_timetree};
use crate::timetree::optimization::relaxed_clock::unit_gammas;
use crate::timetree::params::{TimetreeContext, TimetreeParams};
use crate::timetree::pre_loop::PreLoopState;
use crate::timetree::round::{RoundInputs, RoundState};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_graph::assign_node_names::restrict_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn initialize(
  params: &TimetreeParams,
  context: &TimetreeContext,
  input_names: &BTreeMap<GraphNodeKey, Option<String>>,
  pre_loop: PreLoopState,
  log: &dyn LogSink,
) -> Result<InitialRound, Report> {
  progress_info!(log, "### TreeTime: initialization");
  let PreLoopState {
    graph,
    branch_lengths,
    branch_model,
    clock_fit: ClockFit {
      model: clock_model,
      points: clock_points,
    },
    outliers,
  } = pre_loop;
  let leaf_bad_branches = bad_leaves(&graph, &context.date_constraints, &outliers);
  let state = RoundState {
    names: restrict_node_names(input_names, &graph),
    gammas: unit_gammas(&graph),
    graph,
    branch_model,
    branch_lengths,
    clock_model,
    clock_points,
    clock_branch_lengths: BTreeMap::new(),
    time_inference: TimeInference::default(),
  };
  let inputs = RoundInputs {
    params,
    context,
    leaf_bad_branches: &leaf_bad_branches,
    outliers: &outliers,
  };
  let mode = coalescent_mode(params.coalescent, params.coalescent_opt, params.coalescent_skyline);
  let skyline_params = SkylineParams {
    n_points: params.skyline_n_points,
    stiffness: params.skyline_stiffness,
    n_std: params.coalescent_confidence,
    ..SkylineParams::default()
  };
  let (time_inference, lineages) = if mode.reads_lineages(params.resolve_polytomies) {
    let (time_inference, lineages) = infer_initial_lineages(mode, &skyline_params, &inputs, &state, log)?;
    (time_inference, Some(lineages))
  } else {
    (run_timetree(&state.time_inference_inputs(&inputs), None, log)?, None)
  };

  let state = RoundState {
    time_inference,
    ..state
  };
  Ok(InitialRound {
    state: state.blend_clock_branch_lengths(CLOCK_BRANCH_LENGTH_UNDAMPED, log),
    leaf_bad_branches,
    outliers,
    coalescent: CoalescentSetup {
      mode,
      skyline_params,
      lineages,
    },
  })
}

pub(crate) struct InitialRound {
  pub state: RoundState,
  pub leaf_bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub outliers: BTreeSet<GraphNodeKey>,
  pub coalescent: CoalescentSetup,
}

fn infer_initial_lineages(
  mode: CoalescentMode,
  skyline_params: &SkylineParams,
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<(TimeInference, LineageModel), Report> {
  let run = |prior: Option<&CoalescentPrior>| run_timetree(&state.time_inference_inputs(inputs), prior, log);

  let gaussian_seed = if mode.wants_prior() {
    gaussian_seed(inputs, state, log)?
  } else {
    None
  };
  let (seed_expectation, prior_free) = if let Some(expectation) = gaussian_seed {
    (expectation, None)
  } else {
    let time_inference = run(None)?;
    (
      LineageExpectation::of_inference(&state.graph, &time_inference)?,
      Some(time_inference),
    )
  };
  let timescale = coalescent_timescale(mode, skyline_params, &seed_expectation, log)?;

  let (time_inference, expectation) = if let Some(time_inference) = prior_free.filter(|_| !mode.wants_prior()) {
    (time_inference, seed_expectation)
  } else {
    let prior = CoalescentPrior::new(&seed_expectation, &timescale.distribution)?;
    let time_inference = run(Some(&prior))?;
    let expectation = LineageExpectation::of_inference(&state.graph, &time_inference)?;
    (time_inference, expectation)
  };
  Ok((time_inference, LineageModel { expectation, timescale }))
}

fn gaussian_seed(
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<Option<Arc<LineageExpectation>>, Report> {
  let Some(gaussian) = timetree_gaussian_clock(inputs, state, log) else {
    return Ok(None);
  };
  progress_info!(
    log,
    "Seeding the coalescent prior from the Gaussian node-time marginals"
  );
  let tip_times = given_times(&state.graph, &inputs.context.date_constraints)?;
  let expectation =
    LineageExpectation::from_gaussian(&state.graph, &gaussian.node_times, &tip_times, inputs.leaf_bad_branches)
      .wrap_err("Failed to compute the expected lineage counts of the Gaussian marginals")?;
  Ok(Some(expectation))
}

pub(crate) fn timetree_gaussian_clock(
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Option<GaussianClock> {
  let clock = gaussian_tree(
    &state.graph,
    &state.branch_lengths,
    inputs.context,
    inputs.leaf_bad_branches,
  )
  .and_then(|tree| tree.clock(state.clock_model.clock_rate()));
  match clock {
    Ok(clock) => {
      let fmt_std = |std: Option<f64>| std.map_or_else(|| "undefined".to_owned(), |std| format!("{std:.6e}"));
      progress_info!(
        log,
        "Gaussian clock: rate std {}, root time std {}",
        fmt_std(clock.rate_std),
        fmt_std(clock.root_time_std)
      );
      Some(clock)
    },
    Err(report) => {
      progress_warn!(log, "The Gaussian timetree is not available: {report}");
      None
    },
  }
}

pub(crate) fn gaussian_start_rate(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  context: &TimetreeContext,
  outliers: &BTreeSet<GraphNodeKey>,
  log: &dyn LogSink,
) -> Result<f64, Report> {
  let excluded = bad_leaves(graph, &context.date_constraints, outliers);
  let rate = gaussian_tree(graph, branch_lengths, context, &excluded)
    .and_then(|tree| tree.rate())
    .wrap_err("The root-to-tip regression gives no positive clock rate, and the Gaussian timetree none either")?;
  progress_warn!(
    log,
    "The root-to-tip regression gives no positive clock rate; starting from the rate of the Gaussian timetree, {rate:.6e}"
  );
  Ok(rate)
}

fn gaussian_tree<'g>(
  graph: &'g Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  context: &TimetreeContext,
  excluded: &BTreeMap<GraphNodeKey, bool>,
) -> Result<GaussianTree<'g>, Report> {
  let tip_dates = given_times(graph, &context.date_constraints)?;
  GaussianTree::new(graph, branch_lengths, context.sequence_length, &tip_dates, excluded)
}
