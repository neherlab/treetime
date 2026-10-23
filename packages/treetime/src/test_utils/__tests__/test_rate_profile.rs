#[cfg(test)]
mod tests {
  use crate::test_utils::{FrozenRateState, kingman_neg_log_density, simulated_node_ages};
  use eyre::Report;
  use helpers::{MU, SPACING, small_case};
  use itertools::{Itertools, izip};
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_rate_profile_frozen_prior_is_the_kingman_density_at_the_true_times() -> Result<(), Report> {
    let state = FrozenRateState::new(&small_case()?)?;

    let graph = &state.case.tree.graph;
    let factors = graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        state
          .prior
          .node_contribution(key, node.outbound().len(), state.node_dates[&key])
      })
      .sum::<Result<f64, Report>>()?;

    let expected = kingman_neg_log_density(&simulated_node_ages(&state.case.tree), state.tc)?;
    pretty_assert_abs_diff_eq!(expected, factors - state.prior.log_normalization(), epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_rate_profile_reference_error_bounds_a_four_times_finer_reference() -> Result<(), Report> {
    let state = FrozenRateState::new(&small_case()?)?;
    let rates = [MU / 2.0, MU, MU * 2.0];

    let coarse = state.profile(&rates, true, SPACING)?;
    let fine = state.profile(&rates, true, SPACING / 4.0)?;

    let unbounded = izip!(&coarse, &fine)
      .filter(|(coarse, fine)| (coarse.reference.log_z - fine.reference.log_z).abs() > coarse.reference.error())
      .map(|(coarse, fine)| {
        (
          coarse.rate,
          coarse.reference.log_z,
          fine.reference.log_z,
          coarse.reference.error(),
        )
      })
      .collect_vec();
    assert!(
      unbounded.is_empty(),
      "log Z at spacing h differs from log Z at h/4 by more than the error estimate at (rate, h, h/4, error) {unbounded:?}"
    );
    Ok(())
  }

  // RESEARCH (not for commit): log Z over the rate under the frozen prior, validation layout
  #[test]
  #[ignore = "research"]
  fn research_rate_profile_shape() -> Result<(), Report> {
    use crate::test_utils::{CaseSpec, ModelFamily};
    use rayon::prelude::*;
    use std::fmt::Write as _;
    use treetime_simulation::coalescent::SampleBatch;
    use treetime_simulation::tc_schedule::TcSchedule;
    let sampling = (0..40)
      .map(|i| SampleBatch {
        tau: f64::from(i) * 0.1,
        count: 10,
      })
      .collect_vec();
    let cases = [1000, 10000]
      .into_iter()
      .cartesian_product([3e-4, 1e-3, 3e-3, 1e-2])
      .cartesian_product(0..5_u64)
      .collect_vec();
    let rows = cases
      .par_iter()
      .map(|&((length, mu), seed)| -> Result<String, Report> {
        let spec = CaseSpec {
          seed,
          sampling: sampling.clone(),
          schedule: TcSchedule::constant(10.0)?,
          model: ModelFamily::RandomGtr,
          mu,
          length,
        };
        let state = FrozenRateState::new(&spec)?;
        let mut out = String::new();
        for i in -16..=16 {
          let x = f64::from(i) * 0.1;
          let log_z = state.log_z(mu * x.exp(), true)?.unwrap_or(f64::NAN);
          writeln!(out, "{length}\t{mu}\t{seed}\t{x:.1}\t{log_z}")?;
        }
        eprintln!("done L={length} mu={mu} seed={seed}");
        Ok(out)
      })
      .collect::<Result<Vec<_>, Report>>()?;
    println!("L\tmu\tseed\tx\tlog_z\n{}", rows.concat());
    Ok(())
  }

  mod helpers {
    use crate::test_utils::{CaseSpec, ModelFamily};
    use eyre::Report;
    use treetime_simulation::coalescent::SampleBatch;
    use treetime_simulation::tc_schedule::TcSchedule;

    pub(super) const MU: f64 = 1e-3;
    pub(super) const SPACING: f64 = 1.0 / 32.0;

    pub(super) fn small_case() -> Result<CaseSpec, Report> {
      Ok(CaseSpec {
        seed: 5,
        sampling: vec![SampleBatch { tau: 0.0, count: 3 }, SampleBatch { tau: 0.5, count: 3 }],
        schedule: TcSchedule::constant(1.0)?,
        model: ModelFamily::RandomGtr,
        mu: MU,
        length: 1000,
      })
    }
  }
}
