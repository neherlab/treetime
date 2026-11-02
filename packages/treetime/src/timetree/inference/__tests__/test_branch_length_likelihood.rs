#[cfg(test)]
mod tests {
  use crate::optimize::indel::IndelTerm;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use crate::timetree::inference::runner::{EPS, GRID_POINTS};
  use eyre::Report;
  use rstest::rstest;
  use treetime_distribution::{BoundaryBehavior, Distribution, NegLog};
  use treetime_grid::MaxGridPoints;
  use treetime_utils::array::ndarray::has_uniform_spacing;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::pretty_assert_abs_diff_eq;

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
      IndelTerm::new(0, 0.0),
      0.1,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    pretty_assert_abs_diff_eq!(helpers::eval(&distribution, t), 0.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_indel_rate_only_mode_on_hard_bound() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    pretty_assert_abs_diff_eq!(t_min, 0.0, epsilon = 1e-12);

    let peak_time = distribution.likely_time()?.expect("distribution has a peak");
    pretty_assert_abs_diff_eq!(peak_time, 0.0, epsilon = 1e-12);

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
  #[case::t_0_05( 0.05)]
  #[case::t_0_2(  0.2)]
  #[case::t_0_5(  0.5)]
  #[trace]
  fn test_branch_length_likelihood_indel_rate_only_matches_poisson_shape(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let expected = helpers::INDEL_RATE * t;
    pretty_assert_abs_diff_eq!(helpers::eval(&distribution, t), expected, epsilon = 1e-10);

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
  #[case::just_past_tmax(  1.0)]
  #[case::far_past_tmax(  15.0)]
  #[trace]
  fn test_branch_length_likelihood_right_soft_tail_extrapolates_poisson_decay(#[case] t: f64) -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    assert!(t > t_max, "query {t} must be beyond t_max {t_max}");

    let expected = helpers::INDEL_RATE * t;
    pretty_assert_abs_diff_eq!(distribution.eval(t)?, expected, epsilon = 1e-10);
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
      IndelTerm::new(indel_count, indel_rate),
      5.0,
      one_mutation,
      GRID_POINTS,
      clock_rate * gamma,
      MaxGridPoints::default(),
    )?;

    let t_mle_bl = indel_count.to_f64_exact() / indel_rate;
    let expected_peak_time = t_mle_bl / (clock_rate * gamma);

    let peak_time = distribution.likely_time()?.expect("distribution has a peak");
    pretty_assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-2);
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
      IndelTerm::new(indel_count, indel_rate),
      5.0,
      one_mutation,
      GRID_POINTS,
      clock_rate * gamma,
      MaxGridPoints::default(),
    )?;

    let t_mle_bl = indel_count.to_f64_exact() / indel_rate;
    let expected_peak_time = t_mle_bl / (clock_rate * gamma);

    let peak_time = distribution.likely_time()?.expect("distribution has a peak");
    pretty_assert_abs_diff_eq!(peak_time, expected_peak_time, epsilon = 1e-2);
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
      IndelTerm::new(0, 0.0),
      0.01,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    pretty_assert_abs_diff_eq!(helpers::eval(&with_zero_indels, t), 0.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_grid_extent_scales_with_branch_length() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(0, 0.0),
      0.1,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    pretty_assert_abs_diff_eq!(t_max, 0.5, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_grid_extent_capped_at_max_branch_length() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(0, 0.0),
      2.0,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    pretty_assert_abs_diff_eq!(t_max, 5.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_divergent_boundary_floors_grid_above_zero() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let one_mutation = 1e-3;
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(1, 1.0),
      1.0,
      one_mutation,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    pretty_assert_abs_diff_eq!(t_min, one_mutation * 0.01, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_finite_boundary_grid_starts_at_zero() -> Result<(), Report> {
    let distribution = helpers::build_indel_rate_only_distribution()?;
    let (t_min, _t_max) = distribution.time_bounds().unwrap();
    pretty_assert_abs_diff_eq!(t_min, 0.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_flat_distribution_keeps_pilot_grid() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(0, 0.0),
      0.1,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    let t = distribution.t();
    assert_eq!(GRID_POINTS, t.len());
    assert!(has_uniform_spacing(&t));
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_grid_holds_target_mass_fraction() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let (indel_count, indel_rate) = (1_usize, 10.0);
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(indel_count, indel_rate),
      1.0,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    assert!(
      distribution.t().len() >= GRID_POINTS,
      "stored grid must hold at least GRID_POINTS points, got {}",
      distribution.t().len()
    );

    let (lo, hi) = distribution.time_bounds().unwrap();
    let fraction = helpers::gamma_mass_fraction_inside(indel_count.to_f64_exact(), indel_rate, lo, hi);
    assert!(
      fraction >= 1.0 - 2.0 * EPS,
      "stored grid holds only {fraction} of the mass, below 1 - 2*EPS = {}",
      1.0 - 2.0 * EPS
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::peak_beyond_divergence_span( (5,  100.0,  1e-3,  1e-3))]
  #[case::peak_near_divergence_span(   (2,  1000.0, 2e-3,  1e-4))]
  #[case::peak_inside_divergence_span( (50, 1000.0, 5e-2,  1e-3))]
  #[trace]
  fn test_branch_length_likelihood_grid_holds_peak_and_decaying_tail(
    #[case] (indel_count, indel_rate, current_branch_length, one_mutation): (usize, f64, f64, f64),
  ) -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(indel_count, indel_rate),
      current_branch_length,
      one_mutation,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;
    let k = indel_count.to_f64_exact();

    let Distribution::Function(function) = distribution.as_ref() else {
      panic!("branch-length distribution must be a Function");
    };
    let slope = function.right_extrap().soft_law().expect("a soft Linear tail carries a law").slope;
    assert!(slope > 0.0, "the right tail must decay, got neg-log slope {slope}");

    let peak_time = distribution.likely_time()?.expect("distribution has a peak");
    let spacing = function.dx();
    assert!(
      (peak_time - k / indel_rate).abs() <= spacing,
      "peak {peak_time} must lie within one grid spacing {spacing} of the Gamma mode {}",
      k / indel_rate
    );

    let (lo, hi) = distribution.time_bounds().unwrap();
    let fraction = helpers::gamma_mass_fraction_inside(k, indel_rate, lo, hi);
    assert!(
      fraction >= 1.0 - 2.0 * EPS,
      "stored grid holds only {fraction} of the mass, below 1 - 2*EPS = {}",
      1.0 - 2.0 * EPS
    );

    let survival_quantile = helpers::gamma_survival_quantile(k, indel_rate, EPS);
    assert!(
      hi <= 2.0 * survival_quantile,
      "stored grid ends at {hi}, beyond twice the EPS survival quantile {survival_quantile}"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::extension_mode_with_weak_decay_before_cap( (5,  2.0, 1e-3), 2.5)]
  #[case::extension_mode_beyond_cap(                 (20, 2.0, 1e-3), 5.0)]
  #[case::pilot_mode_with_weak_decay_before_cap(     (5,  2.0, 1.0),  2.5)]
  #[case::pilot_mode_beyond_cap(                     (20, 2.0, 1.0),  5.0)]
  #[trace]
  fn test_branch_length_likelihood_grid_without_decay_before_max_branch_length_ends_hard_at_it(
    #[case] (indel_count, indel_rate, current_branch_length): (usize, f64, f64),
    #[case] expected_peak: f64,
  ) -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(indel_count, indel_rate),
      current_branch_length,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    let Distribution::Function(function) = distribution.as_ref() else {
      panic!("branch-length distribution must be a Function");
    };
    let right = function.right_extrap();
    assert!(
      matches!(right, BoundaryBehavior::Hard),
      "right boundary must be hard at the largest branch length, got {right:?}"
    );

    let (_t_min, t_max) = distribution.time_bounds().unwrap();
    pretty_assert_abs_diff_eq!(t_max, 5.0, epsilon = 1e-12);

    let peak_time = distribution.likely_time()?.expect("distribution has a peak");
    let spacing = function.dx();
    assert!(
      (peak_time - expected_peak).abs() <= spacing,
      "peak {peak_time} must lie within one grid spacing {spacing} of {expected_peak}"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::mode_beyond_cap(              (20, 2.0))]
  #[case::mode_with_weak_decay(         (5,  2.0))]
  #[case::mode_on_hard_bound_slow_decay((0,  1.0))]
  #[trace]
  fn test_branch_length_likelihood_grid_without_decay_before_max_branch_length_ignores_start_length(
    #[case] (indel_count, indel_rate): (usize, f64),
  ) -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let support = |current_branch_length: f64| -> Result<((f64, f64), bool), Report> {
      let distribution = compute_branch_length_distribution(
        &contributions,
        IndelTerm::new(indel_count, indel_rate),
        current_branch_length,
        1e-3,
        GRID_POINTS,
        1.0,
        MaxGridPoints::default(),
      )?;
      let Distribution::Function(function) = distribution.as_ref() else {
        panic!("branch-length distribution must be a Function");
      };
      let hard_right = matches!(function.right_extrap(), BoundaryBehavior::Hard);
      Ok((distribution.time_bounds().unwrap(), hard_right))
    };

    let ((short_min, short_max), short_hard) = support(1e-3)?;
    let ((long_min, long_max), long_hard) = support(1.0)?;

    assert!(short_hard && long_hard, "both right edges must be hard, got {short_hard} and {long_hard}");
    pretty_assert_abs_diff_eq!(short_min, long_min, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(short_max, 5.0, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(long_max, 5.0, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_branch_length_likelihood_constant_likelihood_at_max_branch_length_keeps_soft_edge() -> Result<(), Report> {
    let contributions: Vec<OptimizationContribution> = vec![];
    let distribution = compute_branch_length_distribution(
      &contributions,
      IndelTerm::new(0, 0.0),
      2.0,
      1e-3,
      GRID_POINTS,
      1.0,
      MaxGridPoints::default(),
    )?;

    let Distribution::Function(function) = distribution.as_ref() else {
      panic!("branch-length distribution must be a Function");
    };
    let right = function.right_extrap();
    assert!(
      matches!(right, BoundaryBehavior::Linear(_)),
      "a constant likelihood keeps its soft right edge, got {right:?}"
    );
    Ok(())
  }

  mod helpers {
    use super::*;
    use itertools::izip;
    use std::sync::Arc;
    use treetime_utils::cast::ToF64Exact;

    pub(super) const INDEL_RATE: f64 = 10.0;

    pub(super) fn eval(distribution: &Distribution<NegLog>, t: f64) -> f64 {
      distribution.eval(t).unwrap_or(0.0)
    }

    pub(super) fn gamma_mass_fraction_inside(k: f64, mu: f64, lo: f64, hi: f64) -> f64 {
      let (t, area) = gamma_trapezoid_areas(k, mu);
      let total: f64 = area.iter().sum();
      let inside: f64 = izip!(&t, &area)
        .filter(|(t, _)| **t >= lo && **t <= hi)
        .map(|(_, a)| a)
        .sum();
      inside / total
    }

    pub(super) fn gamma_survival_quantile(k: f64, mu: f64, eps: f64) -> f64 {
      let (t, area) = gamma_trapezoid_areas(k, mu);
      let total: f64 = area.iter().sum();
      let survival = area.iter().rev().scan(0.0, |above, a| {
        *above += a;
        Some(*above)
      });
      let index_from_end = survival.take_while(|above| *above <= eps * total).count();
      t[t.len() - 1 - index_from_end]
    }

    fn gamma_trapezoid_areas(k: f64, mu: f64) -> (Vec<f64>, Vec<f64>) {
      const N: usize = 2_000_001;
      let t_far = (k + 60.0) / mu;
      let dt = t_far / (N.to_f64_exact() - 1.0);
      let density = |t: f64| t.powf(k) * (-mu * t).exp();
      (0..N)
        .map(|i| {
          let t = i.to_f64_exact() * dt;
          let weight = if i == 0 || i == N - 1 { 0.5 } else { 1.0 };
          (t, weight * density(t) * dt)
        })
        .unzip()
    }

    pub(super) fn build_indel_rate_only_distribution() -> Result<Arc<Distribution<NegLog>>, Report> {
      let contributions: Vec<OptimizationContribution> = vec![];
      compute_branch_length_distribution(
        &contributions,
        IndelTerm::new(0, INDEL_RATE),
        1.0,
        1e-3,
        GRID_POINTS,
        1.0,
        MaxGridPoints::default(),
      )
    }
  }
}
