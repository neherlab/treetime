use crate::cancel::Cancel;
use crate::coalescent::coalescent::CoalescentModel;
use crate::progress::{LogSink, StageSink};
use crate::progress_info;
use crate::timetree::coalescent_timescale::{CoalescentSetup, CoalescentTimescale, coalescent_timescale};
use crate::timetree::convergence::metrics::IterationClock;
use crate::timetree::convergence::optimizer::{IterationContext, TimetreeOptimizer, TraceSink};
use crate::timetree::round::{RoundInputs, RoundState, refinement_round};
use eyre::{Report, WrapErr};
use treetime_utils::sync::random::get_random_number_generator;

pub(crate) fn run_refinement_loop(
  inputs: &RoundInputs<'_>,
  coalescent: &CoalescentSetup,
  mut timescale: CoalescentTimescale,
  mut state: RoundState,
  trace_sink: Option<Box<dyn TraceSink + '_>>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(RoundState, CoalescentTimescale), Report> {
  cancel.check()?;
  stages.report("Optimization", 0.3, "");
  progress_info!(log, "### TreeTime: Optimisation rounds");
  let params = inputs.params;
  let mut optimizer = TimetreeOptimizer::new(params.max_iter, false);
  if let Some(sink) = trace_sink {
    optimizer = optimizer.with_trace_sink(sink);
  }
  let max_iter = params.max_iter;

  let seed = params.seed.unwrap_or_else(rand::random);
  if params.resolve_polytomies {
    progress_info!(
      log,
      "Polytomy resolution is stochastic; seed {seed} (pass --seed to reproduce this run)"
    );
  }
  let mut rng = get_random_number_generator(Some(seed));

  while let Some(IterationContext { i }) = optimizer.next_iter(log) {
    cancel.check()?;
    #[expect(
      clippy::as_conversions,
      reason = "iteration counts are far below 2^53, so the conversion to f64 is exact"
    )]
    let iter_fraction = 0.3 + 0.5 * (i as f64 / max_iter as f64);
    stages.report(
      "Optimization",
      iter_fraction,
      &format!("iteration {}/{max_iter}", i + 1),
    );

    if coalescent.mode.is_optimized() {
      timescale = coalescent_timescale(
        coalescent.mode,
        &state.graph,
        &coalescent.skyline_params,
        &state.time_inference.coalescent_node_times()?,
        &state.names,
        log,
      )?;
    }
    let coalescent_model = CoalescentModel::new(&coalescent.lineage_counts, &timescale.distribution)?;
    let merger_rate = coalescent_model.branch_merger_rate_schedule(&timescale.schedule)?;

    let iteration_clock = IterationClock::of(&state.clock_model);
    let (next_state, outcome) = refinement_round(
      inputs,
      &merger_rate,
      coalescent.prior_wanted().then_some(&coalescent_model),
      state,
      &mut rng,
      log,
    )
    .wrap_err_with(|| format!("When running round {i}"))?;
    state = next_state;

    optimizer
      .record(
        outcome.sequence_changes,
        outcome.topology.resolved_nodes(),
        outcome.time_change,
        &state.graph,
        &state.branch_model,
        &state.time_inference,
        coalescent.prior_wanted().then_some(&timescale.distribution),
        iteration_clock,
        &state.names,
        log,
      )
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  Ok((state, timescale))
}
