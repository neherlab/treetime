#[cfg(test)]
mod tests {
  use crate::test_utils::{CaseSpec, recovery_check, require_passed, run_replicates, timetree_cell};
  use eyre::Report;
  use helpers::{
    CLOCK_RATE_BASE_SEED, CLOCK_RATE_REPLICATES, DATE_COVERAGE_BASE_SEED, DATE_COVERAGE_REPLICATES,
    LOW_SIGNAL_BASE_SEED, LOW_SIGNAL_MU, LOW_SIGNAL_REPLICATES, MARGIN, SUBSTITUTIONS_TC_BASE_SEED,
    SUBSTITUTIONS_TC_REPLICATES, TC_BASE_SEED, TC_REPLICATES, clock_rate_error, date_coverage_error,
    estimated_rate_and_tc_params, substitutions_tc_rate_difference, tc_error, true_tc_and_rate_params,
  };
  use treetime_simulation::FALSE_FAIL_RATE;

  #[test]
  fn test_sim_timetree_date_intervals_cover_the_true_dates() -> Result<(), Report> {
    let errors = run_replicates(DATE_COVERAGE_BASE_SEED, DATE_COVERAGE_REPLICATES, |seed| {
      date_coverage_error(&timetree_cell(seed)?, &true_tc_and_rate_params())
    })?;

    require_passed(&[recovery_check(
      "fraction of the true internal node dates inside their 90% interval, minus 0.9, at the true Tc and clock rate",
      &errors,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  #[test]
  #[ignore = "kb/issues/M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md"]
  fn test_sim_timetree_estimated_clock_rate_recovers_the_true_rate() -> Result<(), Report> {
    let errors = run_replicates(CLOCK_RATE_BASE_SEED, CLOCK_RATE_REPLICATES, |seed| {
      clock_rate_error(&timetree_cell(seed)?)
    })?;

    require_passed(&[recovery_check(
      "log ratio of the estimated and the true clock rate, with an estimated constant Tc",
      &errors,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  #[test]
  fn test_sim_timetree_constant_tc_estimate_recovers_the_tc_of_the_true_node_times() -> Result<(), Report> {
    let errors = run_replicates(TC_BASE_SEED, TC_REPLICATES, |seed| tc_error(&timetree_cell(seed)?))?;

    require_passed(&[recovery_check(
      "log ratio of the estimated constant Tc and the maximum likelihood Tc of the true node times (the integral of \
       k(k - 1) / 2 over time, for k lineages, divided by the number of mergers), at the true clock rate",
      &errors,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  #[test]
  fn test_sim_timetree_tc_in_substitutions_per_site_settles_with_the_estimated_rate() -> Result<(), Report> {
    let differences = run_replicates(SUBSTITUTIONS_TC_BASE_SEED, SUBSTITUTIONS_TC_REPLICATES, |seed| {
      substitutions_tc_rate_difference(&timetree_cell(seed)?)
    })?;

    require_passed(&[recovery_check(
      "log ratio of the clock rates estimated on the same case with the true Tc given in substitutions per site at \
       the true rate and with the true Tc given in years; each run with the value in substitutions per site also \
       checks that its Tc in years is that value divided by its final clock rate",
      &differences,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  #[test]
  #[ignore = "kb/issues/M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md"]
  fn test_sim_timetree_date_intervals_cover_the_true_dates_at_low_signal() -> Result<(), Report> {
    let errors = run_replicates(LOW_SIGNAL_BASE_SEED, LOW_SIGNAL_REPLICATES, |seed| {
      let spec = CaseSpec {
        mu: LOW_SIGNAL_MU,
        ..timetree_cell(seed)?
      };
      date_coverage_error(&spec, &estimated_rate_and_tc_params())
    })?;

    require_passed(&[recovery_check(
      "fraction of the true internal node dates inside their 90% interval, minus 0.9, at mu = 3e-4 with an estimated \
       clock rate and Tc; the intervals leave out the rate uncertainty, because confidence is off",
      &errors,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::progress::NoopProgress;
    use crate::test_utils::{
      CaseSpec, SimulatedCase, TIMETREE_CELL_MU, TIMETREE_CELL_TC, match_internal_nodes, pair_exposure,
      simulated_node_ages, timetree_cell_params,
    };
    use crate::timetree::coalescent_timescale::{CoalescentMode, FixedTimescale};
    use crate::timetree::confidence::CI_FRACTION;
    use crate::timetree::params::TimetreeParams;
    use crate::timetree::pipeline::{self, TimetreeOutput};
    use approx::Ulps;
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use treetime_simulation::coalescent::SimulatedTree;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::collections::container::get_exactly_one;
    use treetime_utils::{make_error, make_report};

    pub(super) const MARGIN: f64 = 0.03;
    pub(super) const DATE_COVERAGE_BASE_SEED: u64 = 1;
    pub(super) const DATE_COVERAGE_REPLICATES: usize = 20446;
    pub(super) const CLOCK_RATE_BASE_SEED: u64 = 1_000_001;
    pub(super) const CLOCK_RATE_REPLICATES: usize = 12920;
    pub(super) const TC_BASE_SEED: u64 = 2_000_001;
    pub(super) const TC_REPLICATES: usize = 3813;
    pub(super) const LOW_SIGNAL_BASE_SEED: u64 = 3_000_001;
    pub(super) const LOW_SIGNAL_REPLICATES: usize = 2267;
    pub(super) const LOW_SIGNAL_MU: f64 = 3e-4;
    pub(super) const SUBSTITUTIONS_TC_BASE_SEED: u64 = 4_000_001;
    pub(super) const SUBSTITUTIONS_TC_REPLICATES: usize = 633;

    pub(super) fn true_tc_and_rate_params() -> TimetreeParams {
      TimetreeParams {
        coalescent: CoalescentMode::Fixed(FixedTimescale::Years(TIMETREE_CELL_TC)),
        clock_rate: Some(TIMETREE_CELL_MU),
        ..timetree_cell_params()
      }
    }

    pub(super) fn estimated_rate_and_tc_params() -> TimetreeParams {
      TimetreeParams {
        coalescent: CoalescentMode::Constant,
        ..timetree_cell_params()
      }
    }

    pub(super) fn date_coverage_error(spec: &CaseSpec, params: &TimetreeParams) -> Result<f64, Report> {
      let case = spec.simulate()?;
      let output = run_timetree(&case, spec, params)?;
      let matched = match_internal_nodes(&case.tree, &output.graph, &output.names)?;
      let intervals = output
        .confidence_intervals
        .as_ref()
        .ok_or_else(|| make_report!("The timetree run gave no confidence intervals"))?;
      let intervals = intervals
        .iter()
        .map(|interval| (interval.key, interval))
        .collect::<BTreeMap<_, _>>();
      let covered = matched
        .iter()
        .map(|(&true_node, node)| -> Result<bool, Report> {
          let name = &output.names[*node];
          let interval = intervals
            .get(node)
            .ok_or_else(|| make_report!("Node {name} has no confidence interval"))?;
          let (lower, upper) = (interval.lower, interval.upper);
          if !(lower.is_finite() && upper.is_finite() && lower <= upper) {
            return make_error!("Node {name} has the invalid confidence interval [{lower}, {upper}]");
          }
          let date = case.tree.date(true_node);
          Ok(lower <= date && date <= upper)
        })
        .process_results(|covered| covered.filter(|&covered| covered).count())?;
      Ok(covered.to_f64_exact() / matched.len().to_f64_exact() - CI_FRACTION)
    }

    pub(super) fn clock_rate_error(spec: &CaseSpec) -> Result<f64, Report> {
      let case = spec.simulate()?;
      let output = run_timetree(&case, spec, &estimated_rate_and_tc_params())?;
      let rate = output.clock_model.clock_rate();
      if !(rate.is_finite() && rate > 0.0) {
        return make_error!("The estimated clock rate is {rate}, so its log ratio to the true rate is undefined");
      }
      Ok((rate / spec.mu).ln())
    }

    pub(super) fn tc_error(spec: &CaseSpec) -> Result<f64, Report> {
      let case = spec.simulate()?;
      let params = TimetreeParams {
        clock_rate: Some(spec.mu),
        ..estimated_rate_and_tc_params()
      };
      let output = run_timetree(&case, spec, &params)?;
      let rows = output
        .coalescent
        .as_ref()
        .ok_or_else(|| make_report!("The timetree run gave no coalescent output"))?
        .rows();
      let tc = get_exactly_one(&rows)?.tc_value;
      Ok((tc / constant_tc_of_true_node_times(&case.tree)?).ln())
    }

    pub(super) fn substitutions_tc_rate_difference(spec: &CaseSpec) -> Result<f64, Report> {
      let case = spec.simulate()?;
      let substitutions = TIMETREE_CELL_TC * spec.mu;
      let in_substitutions = run_with_fixed_tc(&case, spec, FixedTimescale::SubstitutionsPerSite(substitutions))?;
      let in_years = run_with_fixed_tc(&case, spec, FixedTimescale::Years(TIMETREE_CELL_TC))?;
      let FixedTcRun { rate, tc } = in_substitutions;
      let settled = substitutions / rate;
      if !Ulps::default().epsilon(0.0).max_ulps(0).eq(&tc, &settled) {
        return make_error!(
          "The Tc in years is {tc}, but {substitutions} substitutions per site at the final clock rate {rate} is {settled}"
        );
      }
      Ok((rate / in_years.rate).ln())
    }

    fn run_with_fixed_tc(
      case: &SimulatedCase,
      spec: &CaseSpec,
      timescale: FixedTimescale,
    ) -> Result<FixedTcRun, Report> {
      let params = TimetreeParams {
        coalescent: CoalescentMode::Fixed(timescale),
        ..timetree_cell_params()
      };
      let output = run_timetree(case, spec, &params)?;
      let rate = output.clock_model.clock_rate();
      if !(rate.is_finite() && rate > 0.0) {
        return make_error!("The estimated clock rate is {rate}, so its log ratio to another rate is undefined");
      }
      let rows = output
        .coalescent
        .as_ref()
        .ok_or_else(|| make_report!("The timetree run gave no coalescent output"))?
        .rows();
      let tc = get_exactly_one(&rows)?.tc_value;
      Ok(FixedTcRun { rate, tc })
    }

    struct FixedTcRun {
      rate: f64,
      tc: f64,
    }

    fn run_timetree(case: &SimulatedCase, spec: &CaseSpec, params: &TimetreeParams) -> Result<TimetreeOutput, Report> {
      Ok(pipeline::run(
        params,
        case.timetree_cell_input(spec.length)?,
        None,
        None,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?)
    }

    fn constant_tc_of_true_node_times(tree: &SimulatedTree) -> Result<f64, Report> {
      let exposure = pair_exposure(&simulated_node_ages(tree))?;
      if exposure.mergers == 0 {
        return make_error!("The true tree has no mergers, so its constant Tc is undefined");
      }
      Ok(exposure.pair_time / exposure.mergers.to_f64_exact())
    }
  }
}
