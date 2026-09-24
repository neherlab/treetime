#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use crate::timetree::inference::runner::{EPS, GRID_POINTS};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_distribution::{BoundaryBehavior, Distribution, NegLog};
  use treetime_utils::array::ndarray::has_uniform_spacing;

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
 0,
 0.0,
 0.1,
 1e-3,
      GRID_POINTS,
 1.0,
 1.0,
    )?;

    assert_abs_diff_eq!(helpers::eval(&distribution, t), 0.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_indel_rate_only_mode_on_hard_bound() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_min, 0.0, epsilon = 1e-12);

    let peak_time = distribution.likely_time().expect("distribution has a peak");
    assert_abs_diff_eq!(peak_time, 0.0, epsilon = 1e-12);

    let Distribution::Function(function) = distribution.as_ref() else {
      panic!("branch-length distribution must be a Function");
    };
    assert!(
      function.y().iter().all(|y| y.is_finite()),
      "no +inf ordinate may be stored on the hard lower bound"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_5(  0.5)]
  #[case::t_2(    2.0)]
  #[case::t_4(    4.0)]
  #[trace]
  fn test_branch_length_likelihood_indel_rate_only_matches_poisson_shape(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let indel_rate = 1.0;
    let expected = indel_rate * t;
    assert_abs_diff_eq!(helpers::eval(&distribution, t), expected, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_right_boundary_is_soft_linear() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let Distribution::Function(function) = distribution.as_ref() else {
      panic!("branch-length distribution must be a Function");
    };
    let right = function.right_extrap();
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

  #[rustfmt::skip]
  #[rstest]
  #[case::just_past_tmax(  9.0)]
  #[case::far_past_tmax(  15.0)]
  #[trace]
  fn test_branch_length_likelihood_right_soft_tail_extrapolates_poisson_decay(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert!(t > t_max, "query {t} must be beyond t_max {t_max}");

    let indel_rate = 1.0;
    let expected = indel_rate * t;
    assert_abs_diff_eq!(distribution.eval(t)?, expected, epsilon = 1e-10);
    Ok(())
  }

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
      5.0,
      one_mutation,
      GRID_POINTS,
      clock_rate,
      gamma,
    )?;

    let t_mle_bl = indel_count as f64 / indel_rate;
    let expected_peak_time = t_mle_bl / (clock_rate * gamma);

    let peak_time = distribution.likely_time().expect("distribution has a peak");
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-2);
    Ok(())
  }

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
      5.0,
      one_mutation,
      GRID_POINTS,
      clock_rate,
      gamma,
    )?;

    let t_mle_bl = indel_count as f64 / indel_rate;
    let expected_peak_time = t_mle_bl / (clock_rate * gamma);

    let peak_time = distribution.likely_time().expect("distribution has a peak");
    assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-2);
    Ok(())
  }

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
 0,
 0.0,
 0.01,
 1e-3,
      GRID_POINTS,
 1.0,
 1.0,
    )?;

    assert_abs_diff_eq!(helpers::eval(&with_zero_indels, t), 0.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_grid_extent_scales_with_branch_length() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(&contributions, 0, 0.0, 0.1, 1e-3, GRID_POINTS, 1.0, 1.0)?;

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_max, 0.5, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_grid_extent_capped_at_max_branch_length() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(&contributions, 0, 0.0, 2.0, 1e-3, GRID_POINTS, 1.0, 1.0)?;

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_max, 5.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_divergent_boundary_floors_grid_above_zero() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let one_mutation = 1e-3;
    let distribution =
      compute_branch_length_distribution(&contributions, 1, 1.0, 1.0, one_mutation, GRID_POINTS, 1.0, 1.0)?;

    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_min, one_mutation * 0.01, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_finite_boundary_grid_starts_at_zero() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;
    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    assert_abs_diff_eq!(t_min, 0.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_flat_distribution_keeps_pilot_grid() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(&contributions, 0, 0.0, 0.1, 1e-3, GRID_POINTS, 1.0, 1.0)?;

    let t = distribution.t();
    assert_eq!(GRID_POINTS, t.len());
    assert!(has_uniform_spacing(&t));
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_grid_holds_target_mass_fraction() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let (indel_count, indel_rate) = (1_usize, 1.0);
    let distribution = compute_branch_length_distribution(
      &contributions,
      indel_count,
      indel_rate,
      1.0,
      1e-3,
      GRID_POINTS,
      1.0,
      1.0,
    )?;

    assert!(
      distribution.t().len() >= GRID_POINTS,
      "stored grid must hold at least GRID_POINTS points, got {}",
      distribution.t().len()
    );

    let (lo, hi) = distribution.time_bounds().unwrap();
    let fraction = helpers::gamma_mass_fraction_inside(indel_count as f64, indel_rate, lo, hi);
    assert!(
      fraction >= 1.0 - 2.0 * EPS,
      "stored grid holds only {fraction} of the mass, below 1 - 2*EPS = {}",
      1.0 - 2.0 * EPS
    );
    Ok(())
  }

  mod helpers {
    use std::sync::Arc;

    use super::*;

    pub(super) fn eval(distribution: &Distribution<NegLog>, t: f64) -> f64 {
      distribution.eval(t).unwrap_or(0.0)
    }

    pub(super) fn gamma_mass_fraction_inside(k: f64, mu: f64, lo: f64, hi: f64) -> f64 {
      const N: usize = 2_000_001;
      let t_far = 60.0;
      let dt = t_far / (N as f64 - 1.0);
      let density = |t: f64| t.powf(k) * (-mu * t).exp();
      let mut total = 0.0;
      let mut inside = 0.0;
      for i in 0..N {
        let t = i as f64 * dt;
        let weight = if i == 0 || i == N - 1 { 0.5 } else { 1.0 };
        let area = weight * density(t) * dt;
        total += area;
        if t >= lo && t <= hi {
          inside += area;
        }
      }
      inside / total
    }

    pub(super) fn build_indel_rate_only_distribution() -> Result<Arc<Distribution<NegLog>>, Report> {
      let contributions: Vec<OptimizationContribution> = vec![];
      compute_branch_length_distribution(&contributions, 0, 1.0, 1.0, 1e-3, GRID_POINTS, 1.0, 1.0)
    }
  }
}
