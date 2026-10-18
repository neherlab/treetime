use crate::clock::clock_regression::ClockFit;
use crate::clock::gaussian::{GaussianClock, gaussian_clock, gaussian_rate};
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::coalescent::skyline::SkylineParams;
use crate::progress::LogSink;
use crate::timetree::coalescent_timescale::{
  CoalescentSetup, CoalescentTimescale, coalescent_mode, coalescent_timescale,
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
use treetime_utils::cast::ToF64Exact;

/// The state that the refinement loop starts from: node times inferred under the coalescent prior, the frozen
/// expected lineage counts and the first coalescent time scale
pub(crate) struct InitialRound {
  pub state: RoundState,
  pub leaf_bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub outliers: BTreeSet<GraphNodeKey>,
  pub coalescent: CoalescentSetup,
  pub timescale: CoalescentTimescale,
}

/// Starts the timetree from the rooted tree and the clock fit of the pre-loop:
///
/// 1. seed the expected lineage counts `kbar(t)` from the Gaussian node-time marginals at the regression rate, or,
///    when the Gaussian timetree degenerates, from a pass of the time inference without a prior;
/// 2. fit the coalescent time scale to the seed;
/// 3. infer the node times under the mean-field prior of the seed and that time scale (with the coalescent prior
///    disabled, the prior-free pass of step 1 is the result);
/// 4. rebuild `kbar(t)` from these posteriors; it stays frozen in the loop.
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
  let run = |prior: Option<&CoalescentPrior>| run_timetree(&state.time_inference_inputs(&inputs), prior, log);
  let expectation_of = |time_inference: &TimeInference| {
    LineageExpectation::from_posterior(&state.graph, &time_inference.posterior, &time_inference.bad_branches)
      .wrap_err("Failed to compute the expected coalescent lineage counts")
  };

  let gaussian_seed = if mode.wants_prior() {
    gaussian_seed(&inputs, &state, log)?
  } else {
    None
  };
  let (seed_expectation, prior_free) = if let Some(expectation) = gaussian_seed {
    (expectation, None)
  } else {
    let time_inference = run(None)?;
    (expectation_of(&time_inference)?, Some(time_inference))
  };
  let timescale = coalescent_timescale(mode, &skyline_params, &seed_expectation, log)?;

  let (time_inference, lineage_expectation) = if let Some(time_inference) = prior_free.filter(|_| !mode.wants_prior()) {
    (time_inference, seed_expectation)
  } else {
    let prior = CoalescentPrior::new(&seed_expectation, &timescale.distribution)?;
    let time_inference = run(Some(&prior))?;
    let expectation = expectation_of(&time_inference)?;
    (time_inference, expectation)
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
      lineage_expectation,
    },
    timescale,
  })
}

/// Expected lineage counts from the Gaussian node-time marginals at the rate of the clock model, or `None` when the
/// Gaussian timetree is not available
fn gaussian_seed(
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<Option<Arc<LineageExpectation>>, Report> {
  let Some(gaussian) = timetree_gaussian_clock(inputs, state, log)? else {
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

/// Gaussian timetree at the rate of the clock model, or `None` when it degenerates
pub(crate) fn timetree_gaussian_clock(
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<Option<GaussianClock>, Report> {
  let tip_dates = given_times(&state.graph, &inputs.context.date_constraints)?;
  match gaussian_clock(
    &state.graph,
    &state.branch_lengths,
    inputs.context.sequence_length.to_f64_exact(),
    &tip_dates,
    inputs.leaf_bad_branches,
    state.clock_model.clock_rate(),
  ) {
    Ok(clock) => {
      progress_info!(
        log,
        "Gaussian clock: rate std {:?}, root time std {:?}",
        clock.rate_std,
        clock.root_time_std
      );
      Ok(Some(clock))
    },
    Err(report) => {
      progress_warn!(log, "The Gaussian timetree is not available: {report}");
      Ok(None)
    },
  }
}

/// Rate of the Gaussian timetree with the root time free, to start from when the root-to-tip regression gives no
/// positive rate; the tips without a date and the `outliers` are left out
pub(crate) fn gaussian_start_rate(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  context: &TimetreeContext,
  outliers: &BTreeSet<GraphNodeKey>,
  log: &dyn LogSink,
) -> Result<f64, Report> {
  let excluded = bad_leaves(graph, &context.date_constraints, outliers);
  let tip_dates = given_times(graph, &context.date_constraints)?;
  let rate = gaussian_rate(
    graph,
    branch_lengths,
    context.sequence_length.to_f64_exact(),
    &tip_dates,
    &excluded,
  )
  .wrap_err("The root-to-tip regression gives no positive clock rate, and the Gaussian timetree none either")?;
  progress_warn!(
    log,
    "The root-to-tip regression gives no positive clock rate; starting from the rate of the Gaussian timetree, {rate:.6e}"
  );
  Ok(rate)
}
