#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::OptimizationContribution;
  use crate::commands::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use treetime_distribution::Distribution;

  const N_GRID: usize = 1000;

  /// Grid spacing in branch-length space for the tests below.
  ///
  /// `create_simple_grid` builds `linspace(one_mutation * 0.1, max_bl, N_GRID)`
  /// with `max_bl = max(center * 3, one_mutation * 10, MAX_BRANCH_TIME * clock_rate)`.
  /// For the test parameters here (`one_mutation = 1e-3`, `clock_rate = 1.0`,
  /// `MAX_BRANCH_TIME = 200.0`), `max_bl = 200` regardless of the center, giving
  /// a uniform spacing of `(200 - 1e-4) / 999 ~ 0.2`. Peak assertions tolerate
  /// one grid cell of slack.
  const GRID_SPACING_BL: f64 = 200.0 / (N_GRID as f64 - 1.0);

  /// Evaluate a distribution at `t` by sampling the underlying function.
  fn eval(distribution: &Distribution, t: f64) -> f64 {
    distribution.eval(t).unwrap_or(0.0)
  }

  /// With no substitution contributions and indel rate zero, every grid point
  /// has log-likelihood zero, so the normalized probability is uniform at 1.0.
  #[test]
  fn test_branch_length_likelihood_no_indels_flat_distribution() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.1,
      /* one_mutation */ 1e-3,
      N_GRID,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    for t in [0.05, 1.0, 10.0, 100.0] {
      assert_abs_diff_eq!(eval(&distribution, t), 1.0, epsilon = 1e-12);
    }
    Ok(())
  }

  /// With no substitution contributions and no observed indels but a positive
  /// indel rate, the Poisson log-likelihood is $-\mu t$, monotonically
  /// decreasing in $t$. The normalized probability in branch-length space is
  /// $\exp(-\mu (t - t_{\min}))$.
  ///
  /// With `clock_rate = gamma = 1.0`, branch-length and time axes coincide, so
  /// interpolated distribution values agree with the closed-form Poisson
  /// expression up to linear interpolation error $\lesssim h^2 |f''| / 8$,
  /// where `h` is the grid spacing.
  #[test]
  fn test_branch_length_likelihood_indel_rate_only_matches_poisson_shape() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let indel_rate = 1.0;
    let one_mutation = 1e-3;
    let clock_rate = 1.0;
    let gamma = 1.0;

    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      indel_rate,
      /* current_branch_length */ 1.0,
      one_mutation,
      N_GRID,
      clock_rate,
      gamma,
    )?;

    // Peak sits at the smallest grid point (log-lh maximized as t -> t_min).
    let t_min = one_mutation * 0.1;
    let expected_peak_time = t_min / (clock_rate * gamma);
    let peak_time = distribution.likely_time().expect("distribution has a peak");
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = GRID_SPACING_BL);

    // Closed-form shape: prob(t) = exp(-mu * (t - t_min)) after normalization.
    // The grid-wise linear interpolation error on exp(-mu t) is bounded by
    // h^2 * mu^2 * exp(-mu t) / 8 = (0.2)^2 / 8 ~ 5e-3 with h = GRID_SPACING_BL.
    // Tolerance 1e-2 is above that bound and orders of magnitude tighter than
    // the "indels ignored" failure mode, where the Poisson term is absent and
    // prob = 1 everywhere.
    //
    // t = 0.01 sits inside the first grid cell [1e-4, 0.2001]; linear
    // interpolation there has O(1%) error, comparable to the tolerance, so
    // we sample only t >= 0.5 where both effects are well separated.
    for t in [0.5, 2.0, 10.0] {
      let expected = (-indel_rate * (t - t_min)).exp();
      assert_abs_diff_eq!(eval(&distribution, t), expected, epsilon = 1e-2);
    }
    Ok(())
  }

  /// With no substitution contributions and `k > 0` observed indels at rate `mu`,
  /// the Poisson log-likelihood peaks at the maximum-likelihood estimate
  /// $\hat{t} = k / \mu$. The distribution in time space peaks at
  /// $\hat{t} / (\text{clock\_rate} \cdot \gamma)$.
  #[test]
  fn test_branch_length_likelihood_indel_mle_peak() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let indel_count: usize = 5;
    let indel_rate = 1.0;
    let clock_rate = 1.0;
    let gamma = 1.0;
    let one_mutation = 1e-3;

    let distribution = compute_branch_length_distribution(
      &contributions,
      indel_count,
      indel_rate,
      /* current_branch_length */ 5.0,
      one_mutation,
      N_GRID,
      clock_rate,
      gamma,
    )?;

    // Poisson MLE in branch-length space: t_mle = k/mu
    let t_mle_bl = indel_count as f64 / indel_rate;
    let expected_peak_time = t_mle_bl / (clock_rate * gamma);

    let peak_time = distribution.likely_time().expect("distribution has a peak");
    // Peak is snapped to the nearest grid point; allow one grid cell of slack.
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = GRID_SPACING_BL);
    Ok(())
  }

  /// Same setup as `test_branch_length_likelihood_indel_mle_peak` with
  /// `gamma = 2.0`: the time-domain peak is compressed by the faster effective
  /// clock rate, confirming the Poisson contribution is added before the
  /// branch-length-to-time conversion.
  #[test]
  fn test_branch_length_likelihood_indel_mle_peak_with_gamma() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let indel_count: usize = 5;
    let indel_rate = 1.0;
    let clock_rate = 1.0;
    let gamma = 2.0;
    let one_mutation = 1e-3;

    let distribution = compute_branch_length_distribution(
      &contributions,
      indel_count,
      indel_rate,
      /* current_branch_length */ 5.0,
      one_mutation,
      N_GRID,
      clock_rate,
      gamma,
    )?;

    let t_mle_bl = indel_count as f64 / indel_rate;
    let expected_peak_time = t_mle_bl / (clock_rate * gamma);

    let peak_time = distribution.likely_time().expect("distribution has a peak");
    let expected_epsilon = GRID_SPACING_BL / (clock_rate * gamma);
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = expected_epsilon);
    Ok(())
  }

  /// Without indels (`indel_count == 0`, `indel_rate == 0`) the Poisson term is
  /// identically zero. The resulting distribution must be identical to one
  /// computed with zero substitutional contributions and no indel arguments,
  /// confirming the new parameters are no-ops on indel-free datasets.
  #[test]
  fn test_branch_length_likelihood_zero_indels_matches_substitution_only() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let with_zero_indels = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.01,
      /* one_mutation */ 1e-3,
      N_GRID,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    for t in [0.001, 0.05, 5.0, 50.0] {
      assert_abs_diff_eq!(eval(&with_zero_indels, t), 1.0, epsilon = 1e-12);
    }
    Ok(())
  }

  /// `compute_branch_length_distribution` rejects non-positive clock rates
  /// before attempting to build the grid; the error message references the
  /// `--clock-rate` flag that the user can set.
  #[test]
  fn test_branch_length_likelihood_rejects_nonpositive_clock_rate() {
    let contributions: Vec<OptimizationContribution> = vec![];
    let result = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.01,
      /* one_mutation */ 1e-3,
      N_GRID,
      /* clock_rate */ -0.5,
      /* gamma */ 1.0,
    );
    let report = result.expect_err("expected negative clock rate to error");
    let rendered = format!("{report:#}");
    assert!(rendered.contains("--clock-rate"), "unexpected error: {rendered}");
  }
}
