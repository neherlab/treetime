use crate::cancel::Cancel;
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::progress::{LogSink, StageSink};
use crate::progress_info;
use crate::timetree::coalescent_timescale::{CoalescentSetup, coalescent_timescale};
use crate::timetree::convergence::metrics::IterationClock;
use crate::timetree::convergence::optimizer::{IterationContext, TimetreeOptimizer, TraceSink};
use crate::timetree::optimization::marginal_clock_rate::RateWindow;
use crate::timetree::round::{RoundInputs, RoundState, refinement_round};
use eyre::{Report, WrapErr};
use std::sync::Arc;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::sync::random::get_random_number_generator;

pub(crate) fn run_refinement_loop(
  inputs: &RoundInputs<'_>,
  mut coalescent: CoalescentSetup,
  mut state: RoundState,
  trace_sink: Option<&mut dyn TraceSink>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(RoundState, CoalescentSetup), Report> {
  cancel.check()?;
  stages.report("Refinement loop", 0.3, "");
  progress_info!(log, "### TreeTime: refinement loop");
  let params = inputs.params;
  let mut optimizer = TimetreeOptimizer::new(params.max_iter);
  if let Some(sink) = trace_sink {
    optimizer = optimizer.with_trace_sink(sink);
  }
  let max_iter = params.max_iter;

  let mut rng = get_random_number_generator(params.seed);
  let mut current = current_expectation(&coalescent, &state)?;

  while let Some(IterationContext { i }) = optimizer.next_iter(log) {
    cancel.check()?;
    let iter_fraction = 0.3 + 0.5 * ((i - 1).to_f64_exact() / max_iter.to_f64_exact());
    stages.report("Refinement loop", iter_fraction, &format!("iteration {i}/{max_iter}"));

    let RoundModels { merger_rate, prior } = round_models(&mut coalescent, current.as_ref(), log)?;

    let iteration_clock = IterationClock::of(&state.clock_model);
    let (next_state, outcome) = refinement_round(
      inputs,
      merger_rate.as_ref(),
      prior.as_ref(),
      if i == 1 {
        RateWindow::FollowSlope
      } else {
        RateWindow::InteriorMaximum
      },
      state,
      &mut rng,
      log,
    )
    .wrap_err_with(|| format!("When running round {i}"))?;
    state = next_state;
    current = current_expectation(&coalescent, &state)?;
    if let (Some(lineages), Some(current), true) = (&mut coalescent.lineages, &current, outcome.topology.changed()) {
      lineages.expectation = Arc::clone(current);
    }
    let expected_log_lh = match (&coalescent.lineages, &current) {
      (Some(lineages), Some(current)) if coalescent.prior_wanted() => {
        Some(lineages.timescale.expected_log_likelihood(current))
      },
      _ => None,
    };

    optimizer
      .record(
        outcome.sequence_changes,
        outcome.topology.resolved_nodes(),
        outcome.time_change,
        &state.graph,
        &state.branch_model,
        &state.time_inference,
        expected_log_lh,
        iteration_clock,
        log,
      )
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  Ok((state, coalescent))
}

fn current_expectation(
  coalescent: &CoalescentSetup,
  state: &RoundState,
) -> Result<Option<Arc<LineageExpectation>>, Report> {
  coalescent
    .lineages
    .as_ref()
    .map(|_| LineageExpectation::of_inference(&state.graph, &state.time_inference))
    .transpose()
}

fn round_models(
  coalescent: &mut CoalescentSetup,
  current: Option<&Arc<LineageExpectation>>,
  log: &dyn LogSink,
) -> Result<RoundModels, Report> {
  let mode = coalescent.mode;
  let prior_wanted = coalescent.prior_wanted();
  let (Some(lineages), Some(current)) = (&mut coalescent.lineages, current) else {
    return Ok(RoundModels {
      merger_rate: None,
      prior: None,
    });
  };
  if mode.is_optimized() {
    lineages.timescale = coalescent_timescale(mode, &coalescent.skyline_params, current, log)?;
  }
  let merger_rate = lineages
    .expectation
    .merger_rate_schedule(&lineages.timescale.distribution)?;
  let prior = prior_wanted
    .then(|| CoalescentPrior::new(&lineages.expectation, &lineages.timescale.distribution))
    .transpose()?;
  Ok(RoundModels {
    merger_rate: Some(merger_rate),
    prior,
  })
}

struct RoundModels {
  merger_rate: Option<PiecewiseConstantFn>,
  prior: Option<CoalescentPrior>,
}
