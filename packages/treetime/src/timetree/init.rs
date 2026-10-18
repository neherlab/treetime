use crate::clock::clock_regression::ClockFit;
use crate::clock::gaussian::{GaussianClock, gaussian_clock};
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
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
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_graph::assign_node_names::restrict_node_names;
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
///    without a sequence length, from a pass of the time inference without a prior;
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
  let (seed_times, seed_expectation, prior_free) = if let Some((times, expectation)) = gaussian_seed {
    (times, expectation, None)
  } else {
    let time_inference = run(None)?;
    let expectation = expectation_of(&time_inference)?;
    (
      time_inference.coalescent_node_times(),
      expectation,
      Some(time_inference),
    )
  };
  let timescale = coalescent_timescale(
    mode,
    &state.graph,
    &skyline_params,
    &seed_times,
    Some(&seed_expectation),
    &state.names,
    log,
  )?;

  let (time_inference, lineage_expectation) = if let Some(time_inference) = prior_free.filter(|_| !mode.wants_prior()) {
    (time_inference, seed_expectation)
  } else {
    let prior = CoalescentPrior::new(&seed_expectation, &timescale.distribution)?;
    let time_inference = run(Some(&prior))?;
    let expectation = expectation_of(&time_inference)?;
    (time_inference, expectation)
  };
  let lineage_counts = compute_lineage_counts(&state.graph, &time_inference.coalescent_node_times())
    .wrap_err("Failed to compute coalescent lineage counts")?;

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
      lineage_counts,
      lineage_expectation,
    },
    timescale,
  })
}

/// Expected lineage counts from the Gaussian node-time marginals at the rate of the clock model, and the point times
/// that lay out the coalescent segments, or `None` when the Gaussian timetree is not available
fn gaussian_seed(
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<Option<(CoalescentNodeTimes, Arc<LineageExpectation>)>, Report> {
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
  let times = point_times(
    &state.graph,
    &expectation.median_times(),
    &tip_times,
    inputs.leaf_bad_branches,
  )?;
  Ok(Some((times, expectation)))
}

/// Gaussian timetree at the rate of the clock model, or `None` when no sequence length is known (no alignment and no
/// `--sequence-length`) or the approximation degenerates
pub(crate) fn timetree_gaussian_clock(
  inputs: &RoundInputs<'_>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<Option<GaussianClock>, Report> {
  let sequence_length = match state.branch_model.sequence_length() {
    0 => inputs.params.sequence_length.unwrap_or(0),
    length => length,
  };
  if sequence_length == 0 {
    debug!("No sequence length: the Gaussian timetree is not computed");
    return Ok(None);
  }
  let tip_dates = given_times(&state.graph, &inputs.context.date_constraints)?;
  match gaussian_clock(
    &state.graph,
    &state.branch_lengths,
    sequence_length.to_f64_exact(),
    &tip_dates,
    inputs.leaf_bad_branches,
    state.clock_model.clock_rate(),
  ) {
    Ok(clock) => {
      progress_info!(
        log,
        "Gaussian clock: rate std {:.4e}, root time std {:.4e}",
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

/// Point times for the lineage counts that lay out the coalescent segments: the dates of the tips and, for an internal
/// node, the median of its seed distribution capped at the times of its children, so that every parent precedes its
/// children; a node without a time is on a bad branch
fn point_times(
  graph: &Graph,
  medians: &BTreeMap<GraphNodeKey, f64>,
  tip_times: &BTreeMap<GraphNodeKey, Option<f64>>,
  leaf_bad_branches: &BTreeMap<GraphNodeKey, bool>,
) -> Result<CoalescentNodeTimes, Report> {
  let mut times: BTreeMap<GraphNodeKey, Option<f64>> = BTreeMap::new();
  graph.iter_depth_first_postorder_forward(|node| {
    let time = if node.is_leaf {
      tip_times.get(&node.key).copied().flatten()
    } else {
      medians.get(&node.key).map(|&median| {
        node
          .child_keys
          .iter()
          .filter_map(|(child, _)| times[child])
          .fold(median, f64::min)
      })
    };
    times.insert(node.key, time);
    Ok(())
  })?;
  Ok(
    times
      .into_iter()
      .map(|(key, time)| {
        let bad_branch = time.is_none() || leaf_bad_branches.get(&key).copied().unwrap_or(false);
        (
          key,
          CoalescentNodeTime {
            time,
            time_dist_likely: time,
            bad_branch,
          },
        )
      })
      .collect(),
  )
}
