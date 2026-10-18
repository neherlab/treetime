use crate::cancel::Cancel;
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::progress::{LogSink, StageSink};
use crate::progress_info;
use crate::timetree::coalescent_timescale::{CoalescentSetup, CoalescentTimescale, coalescent_timescale};
use crate::timetree::convergence::metrics::IterationClock;
use crate::timetree::convergence::optimizer::{IterationContext, TimetreeOptimizer, TraceSink};
use crate::timetree::round::{RoundInputs, RoundState, refinement_round};
use eyre::{Report, WrapErr};
use std::sync::Arc;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::sync::random::get_random_number_generator;

pub(crate) fn run_refinement_loop(
  inputs: &RoundInputs<'_>,
  mut coalescent: CoalescentSetup,
  mut timescale: CoalescentTimescale,
  mut state: RoundState,
  trace_sink: Option<&mut dyn TraceSink>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(RoundState, CoalescentSetup, CoalescentTimescale), Report> {
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
  let expectation_of = |state: &RoundState| {
    LineageExpectation::from_posterior(
      &state.graph,
      &state.time_inference.posterior,
      &state.time_inference.bad_branches,
    )
    .wrap_err("Failed to compute the expected coalescent lineage counts")
  };
  let mut current = expectation_of(&state)?;

  while let Some(IterationContext { i }) = optimizer.next_iter(log) {
    cancel.check()?;
    let iter_fraction = 0.3 + 0.5 * (i.to_f64_exact() / max_iter.to_f64_exact());
    stages.report(
      "Refinement loop",
      iter_fraction,
      &format!("iteration {}/{max_iter}", i + 1),
    );

    if coalescent.mode.is_optimized() {
      timescale = coalescent_timescale(coalescent.mode, &coalescent.skyline_params, &current, log)?;
    }
    let merger_rate = coalescent
      .lineage_expectation
      .merger_rate_schedule(&timescale.distribution)?;
    let coalescent_prior = CoalescentPrior::new(&coalescent.lineage_expectation, &timescale.distribution)?;

    let iteration_clock = IterationClock::of(&state.clock_model);
    let (next_state, outcome) = refinement_round(
      inputs,
      &merger_rate,
      coalescent.prior_wanted().then_some(&coalescent_prior),
      state,
      &mut rng,
      log,
    )
    .wrap_err_with(|| format!("When running round {i}"))?;
    state = next_state;
    current = expectation_of(&state)?;
    if outcome.topology.changed() {
      // the mergers created by polytomy resolution have posteriors now; the expectation takes them in and stays
      // fixed again
      coalescent = CoalescentSetup {
        lineage_expectation: Arc::clone(&current),
        ..coalescent
      };
    }

    optimizer
      .record(
        outcome.sequence_changes,
        outcome.topology.resolved_nodes(),
        outcome.time_change,
        &state.graph,
        &state.branch_model,
        &state.time_inference,
        coalescent
          .prior_wanted()
          .then(|| timescale.expected_log_likelihood(&current)),
        iteration_clock,
        log,
      )
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  Ok((state, coalescent, timescale))
}
