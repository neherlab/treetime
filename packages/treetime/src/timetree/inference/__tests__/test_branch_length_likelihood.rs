#[cfg(test)]
mod tests {
  use crate::partition::optimization_contribution::OptimizationContribution;
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use crate::timetree::inference::runner::BRANCH_GRID_SIZE;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_distribution::{BoundaryBehavior, Distribution, NegLog};
  use treetime_utils::array::ndarray::has_uniform_spacing;

  /// With `current_branch_length = 0.1` the grid spans `[1e-5, 0.5]`; the test
  /// points stay inside that support so the assertion exercises grid values
  /// rather than constant extrapolation beyond the last point.
  ///
  /// A flat likelihood is probability one at every point, whose neg-log ordinate is `-ln(1) = 0`,
  /// so `eval` (which returns the stored ordinate) reads `0.0`, not `1.0`.
  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_001( 0.001)]
  #[case::t_0_05(  0.05)]
  #[case::t_0_2(   0.2)]
  #[case::t_0_49(  0.49)]
  #[trace]
  fn test_branch_length_likelihood_no_indels_flat_distribution(#[case] t: f64) -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.1,
      /* one_mutation */ 1e-3,
      BRANCH_GRID_SIZE,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    assert_abs_diff_eq!(helpers::eval(&distribution, t), 0.0, epsilon = 1e-12);
    Ok(())
  }

  /// With no substitution contributions and no observed indels but a positive
  /// indel rate, the Poisson log-likelihood peaks at $t_{\min}$.
  ///
  /// The grid floor is `min_bl = one_mutation * 0.01`. A pure indel-rate
  /// likelihood decreases monotonically in branch length, so its peak sits at
  /// that floor.
  #[test]
  fn test_branch_length_likelihood_indel_rate_only_peak_at_t_min() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let t_min = 1e-3 * 0.01;
    let expected_peak_time = t_min;
    let peak_time = distribution.likely_time().expect("distribution has a peak");
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-10);

    Ok(())
  }

  /// With no substitution contributions and no observed indels but a positive
  /// indel rate, the normalized probability in branch-length space is
  /// $\exp(-\mu (t - t_{\min}))$, whose neg-log ordinate is $\mu (t - t_{\min})$.
  ///
  /// With `clock_rate = gamma = 1.0`, branch-length and time axes coincide. On the neg-log axis the
  /// stored ordinate is exactly linear in $t$, so grid interpolation reproduces it with no
  /// second-order error: `eval` returns $\mu (t - t_{\min})$ to machine precision.
  ///
  /// The test points stay inside the support so the assertion exercises the
  /// interpolated grid rather than constant extrapolation.
  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_5(  0.5)]
  #[case::t_2(    2.0)]
  #[case::t_4(    4.0)]
  #[trace]
  fn test_branch_length_likelihood_indel_rate_only_matches_poisson_shape(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let indel_rate = 1.0;
    let t_min = 1e-3 * 0.01;
    let expected = indel_rate * (t - t_min);
    assert_abs_diff_eq!(helpers::eval(&distribution, t), expected, epsilon = 1e-10);

    Ok(())
  }

  /// The right boundary is a *soft* `Linear` tail, not the `Error` default: the branch-length
  /// likelihood keeps decaying past the grid's `t_max`, so that edge is where gridding stopped, not
  /// a hard cutoff. A pure indel-rate likelihood decreases rightward, so the fitted neg-log slope is
  /// strictly positive (a genuine, finite-mass decaying tail rather than a flat clamp).
  #[test]
  fn test_branch_length_likelihood_right_boundary_is_soft_linear() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let Distribution::Function(function) = distribution.as_ref() else {
      panic!("branch-length distribution must be a Function");
    };
    let right = function.grid_fn().right_extrap();
    assert!(
      matches!(right, BoundaryBehavior::Linear(_)),
      "right boundary must be a soft Linear tail, got {right:?}"
    );

    let slope = right.soft_law().expect("a soft Linear tail carries a law").slope;
    assert!(
      slope > 0.0,
      "a decaying right tail must have positive neg-log slope, got {slope}"
    );
    Ok(())
  }

  /// Evaluating an inferred branch longer than the grid returns the soft tail's extrapolated neg-log
  /// ordinate instead of failing the `Error` default. Under the `Error` default `eval` would error
  /// and `compute_positional_log_lh` would silently drop the edge; the soft tail keeps it evaluable.
  ///
  /// For the pure indel-rate likelihood the neg-log ordinate is exactly linear, `mu * (t - t_min)`
  /// (the closed-form Poisson decay), so the tail fit reproduces it past `t_max`. The oracle is that
  /// analytic Poisson ordinate, not the system under test.
  #[rustfmt::skip]
  #[rstest]
  #[case::just_past_tmax( 6.0)]
  #[case::far_past_tmax( 10.0)]
  #[trace]
  fn test_branch_length_likelihood_right_soft_tail_extrapolates_poisson_decay(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    // t_max = min(max(1.0 * 5, 1e-3 * 10), 5.0) = 5.0; both query points sit beyond it.
    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert!(t > t_max, "query {t} must be beyond t_max {t_max}");

    let indel_rate = 1.0;
    let t_min = 1e-3 * 0.01;
    let expected = indel_rate * (t - t_min);
    assert_abs_diff_eq!(distribution.eval(t)?, expected, epsilon = 1e-10);
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
      BRANCH_GRID_SIZE,
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
      BRANCH_GRID_SIZE,
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

  /// With `current_branch_length = 0.01` the grid spans `[1e-5, 0.05]`; the
  /// test points stay inside that support so the assertion exercises grid values
  /// rather than constant extrapolation beyond the last point.
  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_001( 0.001)]
  #[case::t_0_01(  0.01)]
  #[case::t_0_03(  0.03)]
  #[case::t_0_05(  0.05)]
  #[trace]
  fn test_branch_length_likelihood_zero_indels_matches_substitution_only(#[case] t: f64) -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let with_zero_indels = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.01,
      /* one_mutation */ 1e-3,
      BRANCH_GRID_SIZE,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    assert_abs_diff_eq!(helpers::eval(&with_zero_indels, t), 0.0, epsilon = 1e-12);
    Ok(())
  }

  /// The grid upper bound follows the ML branch length while the branch stays
  /// below `MAX_BRANCH_LENGTH / 5`: the peak-scaled bound `center * 5` governs
  /// the extent. With `clock_rate = gamma = 1` the time axis equals the
  /// branch-length axis, so the upper time bound is `center * 5`.
  #[test]
  fn test_branch_length_likelihood_grid_extent_scales_with_branch_length() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.1,
      /* one_mutation */ 1e-3,
      BRANCH_GRID_SIZE,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    // max_bl = min(max(0.1 * 5, 1e-3 * 10), 5.0) = min(0.5, 5.0) = 0.5
    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_max, 0.5, epsilon = 1e-12);
    Ok(())
  }

  /// The grid upper bound is capped at `MAX_BRANCH_LENGTH = 5.0` once the
  /// peak-scaled bound `center * 5` would exceed it (branch above 1 sub/site).
  #[test]
  fn test_branch_length_likelihood_grid_extent_capped_at_max_branch_length() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 2.0,
      /* one_mutation */ 1e-3,
      BRANCH_GRID_SIZE,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    // max_bl = min(max(2.0 * 5, 1e-3 * 10), 5.0) = min(10.0, 5.0) = 5.0
    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_max, 5.0, epsilon = 1e-12);
    Ok(())
  }

  /// The grid floor is `one_mutation * 0.01`, independent of branch length. With
  /// `clock_rate = gamma = 1` the lower time bound equals that floor.
  #[test]
  fn test_branch_length_likelihood_grid_lower_bound_from_one_mutation() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let one_mutation = 1e-3;
    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.1,
      one_mutation,
      BRANCH_GRID_SIZE,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_min, one_mutation * 0.01, epsilon = 1e-12);
    Ok(())
  }

  /// The distribution carries exactly the requested number of grid points,
  /// spaced uniformly across `[min_bl, max_bl]`.
  #[test]
  fn test_branch_length_likelihood_grid_has_uniform_requested_resolution() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      /* indel_count */ 0,
      /* indel_rate */ 0.0,
      /* current_branch_length */ 0.1,
      /* one_mutation */ 1e-3,
      BRANCH_GRID_SIZE,
      /* clock_rate */ 1.0,
      /* gamma */ 1.0,
    )?;

    let t = distribution.t();
    assert_eq!(t.len(), BRANCH_GRID_SIZE);
    assert!(has_uniform_spacing(&t));
    Ok(())
  }

  mod helpers {
    use std::sync::Arc;

    use super::*;

    pub fn eval(distribution: &Distribution<NegLog>, t: f64) -> f64 {
      distribution.eval(t).unwrap_or(0.0)
    }

    pub fn build_indel_rate_only_distribution() -> Result<Arc<Distribution<NegLog>>, Report> {
      let contributions: Vec<OptimizationContribution> = vec![];
      compute_branch_length_distribution(
        &contributions,
        /* indel_count */ 0,
        /* indel_rate */ 1.0,
        /* current_branch_length */ 1.0,
        /* one_mutation */ 1e-3,
        BRANCH_GRID_SIZE,
        /* clock_rate */ 1.0,
        /* gamma */ 1.0,
      )
    }
  }
}
