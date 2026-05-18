#[cfg(test)]
mod tests {
  use crate::representation::partition::optimization_contribution::OptimizationContribution;
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_distribution::Distribution;

  const N_GRID: usize = 1000;

  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_05(  0.05)]
  #[case::t_1(     1.0)]
  #[case::t_10(   10.0)]
  #[case::t_100( 100.0)]
  #[trace]
  fn test_branch_length_likelihood_no_indels_flat_distribution(#[case] t: f64) -> Result<(), Report> {
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

    assert_abs_diff_eq!(helpers::eval(&distribution, t), 1.0, epsilon = 1e-12);
    Ok(())
  }

  /// With no substitution contributions and no observed indels but a positive
  /// indel rate, the Poisson log-likelihood peaks at $t_{\min}$.
  #[test]
  fn test_branch_length_likelihood_indel_rate_only_peak_at_t_min() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let t_min = 1e-3 * 0.1;
    let expected_peak_time = t_min;
    let peak_time = distribution.likely_time().expect("distribution has a peak");
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-10);

    Ok(())
  }

  /// With no substitution contributions and no observed indels but a positive
  /// indel rate, the normalized probability in branch-length space is
  /// $\exp(-\mu (t - t_{\min}))$.
  ///
  /// With `clock_rate = gamma = 1.0`, branch-length and time axes coincide, so
  /// interpolated distribution values agree with the closed-form Poisson
  /// expression up to linear interpolation error $\lesssim h^2 |f''| / 8$,
  /// where `h` is the grid spacing (~0.2, giving error bound ~5e-3).
  ///
  /// t >= 0.5 only: the first grid cell [1e-4, 0.2] has O(1%) interpolation
  /// error comparable to the tolerance.
  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_5(  0.5)]
  #[case::t_2(    2.0)]
  #[case::t_10(  10.0)]
  #[trace]
  fn test_branch_length_likelihood_indel_rate_only_matches_poisson_shape(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let indel_rate = 1.0;
    let t_min = 1e-3 * 0.1;
    let expected = (-indel_rate * (t - t_min)).exp();
    assert_abs_diff_eq!(helpers::eval(&distribution, t), expected, epsilon = 1e-2);

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
    // Peak snapped to nearest grid point; measured error 5.1e-3.
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-2);
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
    // Peak snapped to nearest grid point; measured error 2.6e-3.
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-2);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_001( 0.001)]
  #[case::t_0_05(  0.05)]
  #[case::t_5(     5.0)]
  #[case::t_50(   50.0)]
  #[trace]
  fn test_branch_length_likelihood_zero_indels_matches_substitution_only(#[case] t: f64) -> Result<(), Report> {
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

    assert_abs_diff_eq!(helpers::eval(&with_zero_indels, t), 1.0, epsilon = 1e-12);
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

  mod helpers {
    use std::sync::Arc;

    use super::*;

    pub fn eval(distribution: &Distribution, t: f64) -> f64 {
      distribution.eval(t).unwrap_or(0.0)
    }

    pub fn build_indel_rate_only_distribution() -> Result<Arc<Distribution>, Report> {
      let contributions: Vec<OptimizationContribution> = vec![];
      compute_branch_length_distribution(
        &contributions,
        /* indel_count */ 0,
        /* indel_rate */ 1.0,
        /* current_branch_length */ 1.0,
        /* one_mutation */ 1e-3,
        N_GRID,
        /* clock_rate */ 1.0,
        /* gamma */ 1.0,
      )
    }
  }
}
