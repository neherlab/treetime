use crate::clock::clock_regression::{ClockFit, ClockFitResult};
use crate::clock::gaussian::{GaussianClock, GaussianTree};
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::coalescent::skyline::SkylineParams;
use crate::progress::LogSink;
use crate::timetree::coalescent_timescale::{CoalescentMode, CoalescentSetup, LineageModel, coalescent_timescale};
use crate::timetree::inference::bad_branches::bad_leaves;
use crate::timetree::inference::result::{TimeInference, given_times};
use crate::timetree::inference::runner::{CLOCK_BRANCH_LENGTH_UNDAMPED, run_timetree};
use crate::timetree::optimization::marginal_clock_rate::RateWindow;
use crate::timetree::optimization::relaxed_clock::unit_gammas;
use crate::timetree::params::{TimetreeContext, TimetreeParams};
use crate::timetree::pre_loop::PreLoopState;
use crate::timetree::round::{RoundInputs, RoundState};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use itertools::{Itertools, MinMaxResult};
use std::collections::BTreeSet;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_utils::cast::ToF64Exact;

const SIGNAL_RATE_FLOOR_FRACTION: f64 = 0.1;

pub(crate) fn initialize(
  params: &TimetreeParams,
  context: &TimetreeContext,
  pre_loop: PreLoopState,
  log: &dyn LogSink,
) -> Result<InitialRound, Report> {
  progress_info!(log, "### TreeTime: initialization");
  let PreLoopState {
    graph,
    names,
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
    names,
    gammas: unit_gammas(&graph),
    graph,
    branch_model,
    branch_lengths,
    clock_model,
    clock_points,
    rate_window: RateWindow::initial(),
    clock_branch_lengths: EdgeMap::new(),
    time_inference: TimeInference::default(),
  };
  let inputs = RoundInputs {
    params,
    context,
    leaf_bad_branches: &leaf_bad_branches,
    outliers: &outliers,
  };
  let mode = params.coalescent;
  let skyline_params = SkylineParams {
    n_points: params.n_skyline,
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
  pub leaf_bad_branches: NodeMap<bool>,
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
  let timescale = coalescent_timescale(
    mode,
    skyline_params,
    &seed_expectation,
    state.clock_model.clock_rate(),
    log,
  )?;

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

pub(crate) fn clock_fit_or_gaussian_rate(
  fit: ClockFitResult,
  floor: StartRateFloor,
  fit_name: &str,
  tree: &GaussianRateTree<'_>,
  log: &dyn LogSink,
) -> Result<ClockFit, Report> {
  fit.into_clock_fit_or(floor.threshold(), || gaussian_start_rate(tree, floor, fit_name, log))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum StartRateFloor {
  Positive,
  Signal(f64),
}

impl StartRateFloor {
  pub(crate) fn of_fit(fit: &ClockFitResult) -> Self {
    let (dates, divs): (Vec<f64>, Vec<f64>) = fit
      .points()
      .iter()
      .filter(|point| !point.is_outlier)
      .filter_map(|point| point.date.map(|date| (date, point.div)))
      .unzip();
    match dates.iter().copied().minmax() {
      MinMaxResult::MinMax(first, last) if last > first => {
        let mean_div = divs.iter().sum::<f64>() / divs.len().to_f64_exact();
        let floor = SIGNAL_RATE_FLOOR_FRACTION * mean_div / (last - first);
        if floor > 0.0 {
          Self::Signal(floor)
        } else {
          Self::Positive
        }
      },
      MinMaxResult::NoElements | MinMaxResult::OneElement(_) | MinMaxResult::MinMax(..) => Self::Positive,
    }
  }

  fn threshold(self) -> f64 {
    match self {
      Self::Positive => 0.0,
      Self::Signal(floor) => floor,
    }
  }

  fn regression_text(self) -> String {
    match self {
      Self::Positive => "the root-to-tip regression gives no positive clock rate".to_owned(),
      Self::Signal(floor) => format!(
        "the root-to-tip regression gives a clock rate of at most {floor:.6e}, a tenth of the mean root-to-tip divergence per year of sampling"
      ),
    }
  }
}

pub(crate) struct GaussianRateTree<'a> {
  pub graph: &'a Graph,
  pub branch_lengths: &'a EdgeMap<Option<f64>>,
  pub context: &'a TimetreeContext,
  pub outliers: &'a BTreeSet<GraphNodeKey>,
}

fn gaussian_start_rate(
  tree: &GaussianRateTree<'_>,
  floor: StartRateFloor,
  fit_name: &str,
  log: &dyn LogSink,
) -> Result<f64, Report> {
  let regression = floor.regression_text();
  let excluded = bad_leaves(tree.graph, &tree.context.date_constraints, tree.outliers);
  let rate = gaussian_tree(tree.graph, tree.branch_lengths, tree.context, &excluded)
    .and_then(|gaussian| gaussian.rate())
    .wrap_err_with(|| format!("{fit_name}: {regression}, and the Gaussian timetree gives no positive rate either"))?;
  progress_warn!(
    log,
    "{fit_name}: {regression}; using the rate of the Gaussian timetree, {rate:.6e}"
  );
  Ok(rate)
}

fn gaussian_tree<'g>(
  graph: &'g Graph,
  branch_lengths: &EdgeMap<Option<f64>>,
  context: &TimetreeContext,
  excluded: &NodeMap<bool>,
) -> Result<GaussianTree<'g>, Report> {
  let tip_dates = given_times(graph, &context.date_constraints)?;
  GaussianTree::new(graph, branch_lengths, context.sequence_length, &tip_dates, excluded)
}
