use crate::cancel::Cancel;
use crate::coalescent::coalescent::CoalescentModel;
use crate::progress::ProgressSink;
use crate::progress_info;
use crate::timetree::coalescent_timescale::{CoalescentTimescale, coalescent_timescale};
use crate::timetree::convergence::metrics::IterationClock;
use crate::timetree::convergence::optimizer::{IterationContext, TimetreeOptimizer, TraceSink};
use crate::timetree::pipeline::{CoalescentSetup, TimetreeParams};
use crate::timetree::round::{RoundInputs, RoundState, refinement_round};
use eyre::{Report, WrapErr};
use treetime_utils::sync::random::get_random_number_generator;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn run_optimization_loop(
  params: &TimetreeParams,
  inputs: &RoundInputs<'_>,
  coalescent: &CoalescentSetup,
  timescale: CoalescentTimescale,
  state: RoundState,
  trace_sink: Option<Box<dyn TraceSink + '_>>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<(RoundState, CoalescentTimescale), Report> {
  cancel.check()?;
  progress.report("Optimization", 0.3, "");
  progress_info!(progress, "### TreeTime: Optimisation rounds");
  let mut optimizer = TimetreeOptimizer::new(params.max_iter, false);
  if let Some(sink) = trace_sink {
    optimizer = optimizer.with_trace_sink(sink);
  }
  let max_iter = params.max_iter;

  let seed = params.seed.unwrap_or_else(rand::random);
  if params.resolve_polytomies {
    progress_info!(
      progress,
      "Polytomy resolution is stochastic; seed {seed} (pass --seed to reproduce this run)"
    );
  }
  let mut rng = get_random_number_generator(Some(seed));

  let mut state = state;
  let mut timescale = timescale;
  while let Some(IterationContext { i }) = optimizer.next_iter(progress) {
    cancel.check()?;
    let iter_fraction = 0.3 + 0.5 * (i as f64 / max_iter as f64);
    progress.report(
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
        progress,
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
      progress,
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
        progress,
      )
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  Ok((state, timescale))
}
