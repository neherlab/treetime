#[cfg(test)]
mod tests {
  use crate::test_utils::{decrease_check, require_passed, run_replicates, upper_recovery_check};
  use eyre::Report;
  use helpers::{
    MARGIN, RATE_BASE_SEED, RATE_REPLICATE_COUNT, SITES_BASE_SEED, SITES_REPLICATE_COUNT, fixed_rate_log_ratio,
    more_sites_log_ratio,
  };
  use treetime_simulation::FALSE_FAIL_RATE;

  #[test]
  fn test_sim_timetree_more_sites_reduce_the_date_error() -> Result<(), Report> {
    let log_ratios = run_replicates(SITES_BASE_SEED, SITES_REPLICATE_COUNT, more_sites_log_ratio)?;

    require_passed(&[decrease_check(
      "log ratio of the internal-node date RMSE at 4000 and at 500 sites",
      &log_ratios,
      FALSE_FAIL_RATE,
    )?])
  }

  #[test]
  fn test_sim_timetree_fixed_true_rate_does_not_worsen_the_dates() -> Result<(), Report> {
    let errors = run_replicates(RATE_BASE_SEED, RATE_REPLICATE_COUNT, fixed_rate_log_ratio)?;

    require_passed(&[upper_recovery_check(
      "log ratio of the internal-node date RMSE with the true rate fixed and with an estimated rate",
      &errors,
      MARGIN,
      FALSE_FAIL_RATE,
    )?])
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::progress::NoopProgress;
    use crate::test_utils::{
      CaseSpec, SimulatedCase, TIMETREE_CELL_LENGTH, TIMETREE_CELL_MU, TIMETREE_CELL_TC, match_internal_nodes,
      timetree_cell, timetree_cell_params,
    };
    use crate::timetree::params::TimetreeParams;
    use crate::timetree::pipeline;
    use eyre::Report;
    use ndarray::Array1;
    use treetime_utils::{make_internal_report, make_report};

    pub(super) const SITES_BASE_SEED: u64 = 1;
    pub(super) const SITES_REPLICATE_COUNT: usize = 15;
    pub(super) const RATE_BASE_SEED: u64 = 100_001;
    pub(super) const RATE_REPLICATE_COUNT: usize = 812;
    pub(super) const MARGIN: f64 = 0.03;
    const FEW_SITES: usize = 500;
    const ALL_SITES: usize = 4000;

    pub(super) fn more_sites_log_ratio(seed: u64) -> Result<f64, Report> {
      let case = CaseSpec {
        length: ALL_SITES,
        ..timetree_cell(seed)?
      }
      .simulate()?;
      let few = date_rmse(&case, &true_cell_params(Some(TIMETREE_CELL_MU)), FEW_SITES)?;
      let all = date_rmse(&case, &true_cell_params(Some(TIMETREE_CELL_MU)), ALL_SITES)?;
      Ok((all / few).ln())
    }

    pub(super) fn fixed_rate_log_ratio(seed: u64) -> Result<f64, Report> {
      let case = timetree_cell(seed)?.simulate()?;
      let fixed = date_rmse(&case, &true_cell_params(Some(TIMETREE_CELL_MU)), TIMETREE_CELL_LENGTH)?;
      let estimated = date_rmse(&case, &true_cell_params(None), TIMETREE_CELL_LENGTH)?;
      Ok((fixed / estimated).ln())
    }

    fn true_cell_params(clock_rate: Option<f64>) -> TimetreeParams {
      TimetreeParams {
        clock_rate,
        coalescent: Some(TIMETREE_CELL_TC),
        ..timetree_cell_params()
      }
    }

    fn date_rmse(case: &SimulatedCase, params: &TimetreeParams, columns: usize) -> Result<f64, Report> {
      let output = pipeline::run(
        params,
        case.timetree_cell_input(columns)?,
        None,
        None,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let matched = match_internal_nodes(&case.tree, &output.graph, &output.names)?;
      let squared_errors = matched
        .iter()
        .map(|(&true_node, &node)| {
          let date = output.node_dates[&node]
            .ok_or_else(|| make_report!("Internal node {} has no date", &output.names[node]))?;
          Ok((date - case.tree.date(true_node)).powi(2))
        })
        .collect::<Result<Array1<f64>, Report>>()?;
      let mean = squared_errors
        .mean()
        .ok_or_else(|| make_internal_report!("The true tree has no internal node"))?;
      Ok(mean.sqrt())
    }
  }
}
