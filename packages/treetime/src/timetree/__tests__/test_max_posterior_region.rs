#[cfg(test)]
mod tests {
  use crate::timetree::max_posterior_region::{
    MaxPosteriorRegion, PeakLocation, PosteriorRegion, max_posterior_region,
  };
  use eyre::Report;
  use helpers::{assert_second_order, function, grid, grid_region, model_density, region_errors, trapezoid};
  use itertools::chain;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use statrs::distribution::{Continuous, ContinuousCDF, Gamma, Normal};
  use std::f64::consts::LN_10;
  use std::iter::{once, repeat_n};
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, MassProfile, NegLog, TailMass};
  use treetime_grid::{HardApproachLaw, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq, prop_assert_abs_diff_eq};

  const HARD: BoundaryBehavior = BoundaryBehavior::Hard;
  const CELLS: [usize; 4] = [100, 200, 400, 800];
  const Z_95: f64 = 1.644_853_626_951_472_2;

  #[test]
  fn test_max_posterior_region_gaussian_converges_to_symmetric_hpd() -> Result<(), Report> {
    let (mean, std_dev) = (2015.0, 0.5);
    let gaussian = |cells: usize| {
      function(
        (mean - 8.0 * std_dev, mean + 8.0 * std_dev),
        cells + 1,
        HARD,
        HARD,
        |x| 3.7 + (x - mean).powi(2) / (2.0 * std_dev * std_dev),
      )
    };
    let expected = (mean - Z_95 * std_dev, mean + Z_95 * std_dev);

    assert_second_order(&region_errors(&CELLS, gaussian, expected)?);
    let region = grid_region(&gaussian(160_000)?, 0.9)?;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(expected.0, region.lower, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(expected.1, region.upper, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_skewed_gamma_has_equal_density_bounds() -> Result<(), Report> {
    let (shape, rate) = (3.0, 2.0);
    let distribution = function((1e-4, 20.0), 200_000, HARD, HARD, |x| {
      -(shape - 1.0) * x.ln() + rate * x
    })?;

    let region = grid_region(&distribution, 0.9)?;

    let gamma = Gamma::new(shape, rate)?;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(0.9, gamma.cdf(region.upper) - gamma.cdf(region.lower), epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(gamma.pdf(region.lower), gamma.pdf(region.upper), epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::level_below_the_second_mode_keeps_main_mode(0.6,  false)]
  #[case::level_between_second_mode_and_dip(           0.8,  false)]
  #[case::level_above_the_dip_adds_second_mode(       0.99, true)]
  #[trace]
  fn test_max_posterior_region_bimodal_region_is_level_set_connected_to_peak(
    #[case] fraction: f64,
    #[case] contains_second_mode: bool,
  ) -> Result<(), Report> {
    let main = Normal::new(5.0, 0.5)?;
    let side_mode = 3.2;
    let side = Normal::new(side_mode, 0.3)?;
    let (main_weight, side_weight) = (0.85, 0.15);
    let neglog = |x: f64| -(main_weight * main.pdf(x) + side_weight * side.pdf(x)).ln();
    let t = grid(-2.0, 10.0, 120_001);
    let y = t.mapv(neglog);
    let distribution = function((-2.0, 10.0), 120_001, HARD, HARD, neglog)?;

    let region = grid_region(&distribution, fraction)?;

    let mixture_cdf = |x: f64| main_weight * main.cdf(x) + side_weight * side.cdf(x);
    let level = model_density(&t, &y, region.upper);
    let lowest_inside = t
      .iter()
      .zip(&y)
      .filter(|(x, _)| **x > region.lower && **x < region.upper)
      .map(|(_, y)| (-y).exp())
      .fold(f64::INFINITY, f64::min);
    assert_eq!(PeakLocation::Interior, region.peak);
    assert_eq!(contains_second_mode, region.lower < side_mode);
    pretty_assert_abs_diff_eq!(fraction, mixture_cdf(region.upper) - mixture_cdf(region.lower), epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(level, model_density(&t, &y, region.lower), epsilon = 1e-15);
    assert!(lowest_inside >= level, "every point inside {region:?} lies above the density level");
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_peak_at_earliest_date_extends_later() -> Result<(), Report> {
    let exponential = |cells: usize| function((2000.0, 2010.0), cells + 1, HARD, HARD, |x| x - 2000.0);
    let truncated_exponential_quantile = 2000.0 - (1.0 - 0.9 * (1.0 - (-10.0_f64).exp())).ln();

    assert_second_order(&region_errors(
      &CELLS,
      exponential,
      (2000.0, truncated_exponential_quantile),
    )?);
    let region = grid_region(&exponential(100_000)?, 0.9)?;
    assert_eq!(PeakLocation::Earliest, region.peak);
    pretty_assert_ulps_eq!(2000.0, region.lower, max_ulps = 0);
    pretty_assert_abs_diff_eq!(truncated_exponential_quantile, region.upper, epsilon = 1e-13);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_peak_at_latest_date_extends_earlier() -> Result<(), Report> {
    let exponential = |cells: usize| function((2000.0, 2010.0), cells + 1, HARD, HARD, |x| 2010.0 - x);
    let truncated_exponential_quantile = 2010.0 + (1.0 - 0.9 * (1.0 - (-10.0_f64).exp())).ln();

    assert_second_order(&region_errors(
      &CELLS,
      exponential,
      (truncated_exponential_quantile, 2010.0),
    )?);
    let region = grid_region(&exponential(100_000)?, 0.9)?;
    assert_eq!(PeakLocation::Latest, region.peak);
    pretty_assert_abs_diff_eq!(truncated_exponential_quantile, region.lower, epsilon = 1e-13);
    pretty_assert_ulps_eq!(2010.0, region.upper, max_ulps = 0);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::grid_within_one_decay_length(   1.0)]
  #[case::grid_within_three_decay_lengths(3.0)]
  #[trace]
  fn test_max_posterior_region_laplace_with_exact_tails_converges_to_log_ten_over_rate(
    #[case] half_width_in_decay_lengths: f64,
  ) -> Result<(), Report> {
    let (t0, rate) = (2015.0, 2.0);
    let half_width = half_width_in_decay_lengths / rate;
    let laplace = |cells: usize| {
      function(
        (t0 - half_width, t0 + half_width),
        cells + 1,
        BoundaryBehavior::Linear(SoftTailLaw { slope: -rate }),
        BoundaryBehavior::Linear(SoftTailLaw { slope: rate }),
        |x| rate * (x - t0).abs(),
      )
    };
    let expected = (t0 - LN_10 / rate, t0 + LN_10 / rate);

    assert_second_order(&region_errors(&CELLS, laplace, expected)?);
    let region = grid_region(&laplace(100_000)?, 0.9)?;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(expected.0, region.lower, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(expected.1, region.upper, epsilon = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::narrow_at_year_2000(  2000.0,  1e6 )]
  #[case::narrow_at_year_0(     0.0,     1e6 )]
  #[case::narrow_at_year_1(     1.0,     1e6 )]
  #[case::narrow_before_year_0( -2000.0, 1e6 )]
  #[case::narrower_at_year_2000(2000.0,  1e9 )]
  #[case::unit_rate_at_year_2000(2000.0, 1.0 )]
  #[case::wide_at_year_2000(    2000.0,  1e-3)]
  #[trace]
  fn test_max_posterior_region_three_point_laplace_is_exact_under_translation_and_scale(
    #[case] t0: f64,
    #[case] rate: f64,
  ) -> Result<(), Report> {
    let half_width = 1.0 / rate;
    let distribution = function(
      (t0 - half_width, t0 + half_width),
      3,
      BoundaryBehavior::Linear(SoftTailLaw { slope: -rate }),
      BoundaryBehavior::Linear(SoftTailLaw { slope: rate }),
      |x| rate * (x - t0).abs(),
    )?;

    let region = grid_region(&distribution, 0.9)?;

    let total_times_rate = 2.0 * trapezoid(1.0, 0.0, 1.0) + 2.0 * (-1.0_f64).exp();
    let half_region = -(0.05 * total_times_rate).ln() / rate;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(t0 - half_region, region.lower, epsilon = 1e-11);
    pretty_assert_abs_diff_eq!(t0 + half_region, region.upper, epsilon = 1e-11);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_laplace_on_one_sided_grid_continues_into_earlier_tail() -> Result<(), Report> {
    let (t0, rate) = (2015.0, 2.0);
    let laplace = |cells: usize| {
      function(
        (t0, t0 + 1.0 / rate),
        cells + 1,
        BoundaryBehavior::Linear(SoftTailLaw { slope: -rate }),
        BoundaryBehavior::Linear(SoftTailLaw { slope: rate }),
        |x| rate * (x - t0),
      )
    };
    let expected = (t0 - LN_10 / rate, t0 + LN_10 / rate);

    assert_second_order(&region_errors(&CELLS, laplace, expected)?);
    let region = grid_region(&laplace(100_000)?, 0.9)?;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(expected.0, region.lower, epsilon = 1e-11);
    pretty_assert_abs_diff_eq!(expected.1, region.upper, epsilon = 1e-11);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_exponential_from_hard_edge_ends_in_soft_tail() -> Result<(), Report> {
    let (t0, rate) = (2015.0, 2.0);
    let exponential = |cells: usize| {
      function(
        (t0, t0 + 1.0 / rate),
        cells + 1,
        HARD,
        BoundaryBehavior::Linear(SoftTailLaw { slope: rate }),
        |x| rate * (x - t0),
      )
    };

    assert_second_order(&region_errors(&CELLS, exponential, (t0, t0 + LN_10 / rate))?);
    let region = grid_region(&exponential(100_000)?, 0.9)?;
    assert_eq!(PeakLocation::Earliest, region.peak);
    pretty_assert_ulps_eq!(t0, region.lower, max_ulps = 0);
    pretty_assert_abs_diff_eq!(t0 + LN_10 / rate, region.upper, epsilon = 1e-11);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_power_law_tails_hold_the_fraction() -> Result<(), Report> {
    let (half_width, b) = (0.1_f64, 2.0_f64);
    let y_edge = -b * (1.0 - half_width).ln();
    let distribution = function(
      (-half_width, half_width),
      3,
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: -1.0, b }),
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 1.0, b }),
      |x| if x == 0.0 { 0.0 } else { y_edge },
    )?;

    let region = grid_region(&distribution, 0.9)?;

    let tail_mass_beyond = |distance_to_hard_edge: f64| distance_to_hard_edge.powf(b + 1.0) / (b + 1.0);
    let total = 2.0 * (trapezoid(half_width, 0.0, y_edge) + tail_mass_beyond(1.0 - half_width));
    let half_region = 1.0 - (0.05 * (b + 1.0) * total).powf((b + 1.0).recip());
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(-half_region, region.lower, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(half_region, region.upper, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_inside_hard_cells_inverts_their_power_law() -> Result<(), Report> {
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values((-1.0, 1.0), array![f64::INFINITY, 0.0, f64::INFINITY])?
        .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: -1.0, b: 1.0 }))?
        .with_right_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 1.0, b: 1.0 }))?,
    );

    let region = grid_region(&distribution, 0.9)?;

    let triangle_cdf = |x: f64| {
      if x < 0.0 {
        (1.0 + x).powi(2) / 2.0
      } else {
        1.0 - (1.0 - x).powi(2) / 2.0
      }
    };
    let half_region = 1.0 - 0.1_f64.sqrt();
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(-half_region, region.lower, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(half_region, region.upper, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(
      0.9,
      triangle_cdf(region.upper) - triangle_cdf(region.lower),
      epsilon = 1e-15
    );
    pretty_assert_abs_diff_eq!(1.0 - region.lower.abs(), 1.0 - region.upper.abs(), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_jump_at_plateau_level_ends_inside_the_plateau() -> Result<(), Report> {
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(
        (0.0, 7.0),
        array![1.0, 0.0, 5.0, 0.5, 0.5, 0.5, 0.5, 30.0],
      )?
      .with_left_extrap(HARD)?
      .with_right_extrap(HARD)?,
    );

    let region = grid_region(&distribution, 0.9)?;

    let (rise, fall, dip, plateau, end) = (
      trapezoid(1.0, 1.0, 0.0),
      trapezoid(1.0, 0.0, 5.0),
      trapezoid(1.0, 5.0, 0.5),
      3.0 * trapezoid(1.0, 0.5, 0.5),
      trapezoid(1.0, 0.5, 30.0),
    );
    let total = rise + fall + dip + plateau + end;
    assert!((rise + fall) / total < 0.9);
    assert!(rise + fall + dip < 0.9 * total && 0.9 * total < rise + fall + dip + plateau);
    let upper_in_plateau = 3.0 + (0.9 * total - rise - fall - dip) / (-0.5_f64).exp();
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_ulps_eq!(0.0, region.lower, max_ulps = 0);
    pretty_assert_abs_diff_eq!(upper_in_plateau, region.upper, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_zero_density_gap_is_crossed_inside_the_next_component() -> Result<(), Report> {
    let function = DistributionFunction::<f64, NegLog>::from_range_values(
      (0.0, 6.0),
      array![1.0, 0.0, 1.0, f64::INFINITY, 0.1, 0.1, 1.0],
    )?
    .with_left_extrap(HARD)?
    .with_right_extrap(HARD)?;

    let region = grid_region(&Distribution::Function(function.clone()), 0.9)?;

    let profile = MassProfile::new(&function)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let total = fractions.total();
    let up_to_flat_cell = fractions.mass_between(0.0, 4.0)?;
    let flat_cell = trapezoid(1.0, 0.1, 0.1);
    assert!(up_to_flat_cell / total < 0.9 && (up_to_flat_cell + flat_cell) / total < 0.9);
    let in_falling_cell = 0.9 * total - up_to_flat_cell - flat_cell;
    let (p_start, p_end) = ((-0.1_f64).exp(), (-1.0_f64).exp());
    let slope = p_end - p_start;
    let upper_in_falling_cell =
      5.0 + 2.0 * in_falling_cell / (p_start + (p_start * p_start + 2.0 * slope * in_falling_cell).sqrt());
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_ulps_eq!(0.0, region.lower, max_ulps = 0);
    pretty_assert_abs_diff_eq!(upper_in_falling_cell, region.upper, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(
      0.9,
      fractions.mass_between(region.lower, region.upper)? / total,
      epsilon = 1e-15
    );
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_flat_hard_approach_tail_is_crossed_by_the_jump_rule() -> Result<(), Report> {
    let t_hard = -10.0;
    let distribution = function(
      (0.0, 2.0),
      3,
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard, b: 0.0 }),
      HARD,
      |x| x,
    )?;

    let region = grid_region(&distribution, 0.9)?;

    let flat_tail = -t_hard;
    let total = flat_tail + trapezoid(1.0, 0.0, 1.0) + trapezoid(1.0, 1.0, 2.0);
    assert!(0.9 * total < flat_tail);
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(-0.9 * total, region.lower, epsilon = 1e-15);
    pretty_assert_ulps_eq!(0.0, region.upper, max_ulps = 0);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_flat_top_holding_the_fraction_extends_from_latest_peak() -> Result<(), Report> {
    let y: Array1<f64> = chain!(once(20.0), repeat_n(0.0, 21), once(20.0)).collect();
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values((0.0, 22.0), y)?
        .with_left_extrap(HARD)?
        .with_right_extrap(HARD)?,
    );

    let region = grid_region(&distribution, 0.9)?;

    let total = 20.0 + 2.0 * trapezoid(1.0, 0.0, 20.0);
    assert!(20.0 / total > 0.9);
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(21.0 - 0.9 * total, region.lower, epsilon = 1e-15);
    pretty_assert_ulps_eq!(21.0, region.upper, max_ulps = 0);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_single_finite_ordinate_spreads_over_its_flat_hard_cell() -> Result<(), Report> {
    let distribution = function(
      (0.0, 1.0),
      2,
      HARD,
      HARD,
      |x| if x == 0.0 { 0.0 } else { f64::INFINITY },
    )?;

    let region = grid_region(&distribution, 0.9)?;

    let flat_cell_quantile = PosteriorRegion {
      lower: 0.0,
      upper: 0.9,
      peak: PeakLocation::Earliest,
    };
    assert_eq!(flat_cell_quantile, region);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_flat_range_resolves_tie_to_latest_date() -> Result<(), Report> {
    let distribution = Distribution::<NegLog>::range((2019.0, 2021.0), 4.0);

    let region = max_posterior_region(&distribution, 0.9)?;

    let expected = MaxPosteriorRegion::Grid(PosteriorRegion {
      lower: 2019.2,
      upper: 2021.0,
      peak: PeakLocation::Latest,
    });
    assert_eq!(Some(expected), region);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_point_distribution_is_a_point_region() -> Result<(), Report> {
    let distribution = Distribution::<NegLog>::point(2020.5, 1.0);

    let region = max_posterior_region(&distribution, 0.9)?;

    let expected = MaxPosteriorRegion::Point;
    assert_eq!(Some(expected), region);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_zero_probability_point_has_no_region() -> Result<(), Report> {
    let distribution = Distribution::<NegLog>::point(2020.5, f64::INFINITY);
    assert_eq!(None, max_posterior_region(&distribution, 0.9)?);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_empty_has_no_region() -> Result<(), Report> {
    assert_eq!(None, max_posterior_region(&Distribution::<NegLog>::empty(), 0.9)?);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_function_without_finite_ordinate_has_no_region() -> Result<(), Report> {
    let distribution = function((0.0, 1.0), 2, HARD, HARD, |_| f64::INFINITY)?;
    assert_eq!(None, max_posterior_region(&distribution, 0.9)?);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_far_zero_density_cells_leave_the_region_unchanged() -> Result<(), Report> {
    let (mean, std_dev) = (2015.0, 0.5);
    let neglog = |x: f64| (x - mean).powi(2) / (2.0 * std_dev * std_dev);
    let range = (mean - 8.0 * std_dev, mean + 8.0 * std_dev);
    let distribution = function(range, 16001, HARD, HARD, neglog)?;
    let distribution_hard_tails = function(range, 16001, HARD, HARD, |x| {
      let y = neglog(x);
      if y > 30.0 { f64::INFINITY } else { y }
    })?;

    let region = grid_region(&distribution, 0.9)?;
    let region_hard_tails = grid_region(&distribution_hard_tails, 0.9)?;

    assert_eq!(PeakLocation::Interior, region_hard_tails.peak);
    pretty_assert_ulps_eq!(region.lower, region_hard_tails.lower, max_ulps = 0);
    pretty_assert_ulps_eq!(region.upper, region_hard_tails.upper, max_ulps = 0);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_rejects_nan_ordinates() -> Result<(), Report> {
    let distribution = Distribution::Function(DistributionFunction::<f64, NegLog>::from_range_values(
      (0.0, 2.0),
      array![0.0, f64::NAN, 1.0],
    )?);
    let result = max_posterior_region(&distribution, 0.9);
    assert_error!(
      result,
      "Cannot normalize a distribution on [0, 2]: its negative log-likelihood values contain NaN"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(    0.0)]
  #[case::one(     1.0)]
  #[case::negative(-0.1)]
  #[case::nan(     f64::NAN)]
  #[trace]
  fn test_max_posterior_region_rejects_fraction_outside_unit_interval(#[case] fraction: f64) {
    let result = max_posterior_region(&Distribution::<NegLog>::range((0.0, 1.0), 0.0), fraction);
    assert_error!(result, format!("Posterior mass fraction must be in (0, 1), got {fraction}"));
  }

  proptest! {
    #[test]
    fn test_prop_max_posterior_region_interior_bounds_have_equal_model_density(
      center in 0.4_f64..0.6,
      width_lower in 0.01_f64..0.1,
      width_upper in 0.01_f64..0.1,
      power in 1.0_f64..3.0,
    ) {
      let t = grid(0.0, 1.0, 2001);
      let shape = |x: f64| {
        let width = if x < center { width_lower } else { width_upper };
        ((x - center).abs() / width).powf(power)
      };
      let y = t.mapv(shape);
      let distribution = function((0.0, 1.0), 2001, HARD, HARD, shape).unwrap();

      let region = grid_region(&distribution, 0.9).unwrap();

      prop_assume!(region.peak == PeakLocation::Interior);
      prop_assert_abs_diff_eq!(model_density(&t, &y, region.lower), model_density(&t, &y, region.upper), epsilon = 1e-13);
    }

    #[test]
    fn test_prop_max_posterior_region_affine_time_change_preserves_region(
      center in 0.2_f64..0.8,
      width in 0.02_f64..0.3,
      scale in 0.1_f64..10.0,
      shift in -3000.0_f64..3000.0,
    ) {
      let y = grid(0.0, 1.0, 1001).mapv(|x| ((x - center) / width).powi(2) + (x - center));
      let distribution = |range: (f64, f64)| {
        Distribution::Function(
          DistributionFunction::<f64, NegLog>::from_range_values(range, y.clone())
            .and_then(|f| f.with_left_extrap(HARD))
            .and_then(|f| f.with_right_extrap(HARD))
            .unwrap(),
        )
      };

      let region = grid_region(&distribution((0.0, 1.0)), 0.9).unwrap();
      let moved = grid_region(&distribution((shift, scale + shift)), 0.9).unwrap();

      prop_assert_eq!(region.peak, moved.peak);
      prop_assert_abs_diff_eq!(scale * region.lower + shift, moved.lower, epsilon = 1e-12);
      prop_assert_abs_diff_eq!(scale * region.upper + shift, moved.upper, epsilon = 1e-12);
    }
  }

  mod helpers {
    use crate::timetree::max_posterior_region::{MaxPosteriorRegion, PosteriorRegion, max_posterior_region};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::Array1;
    use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
    use treetime_utils::make_error;

    pub(super) fn grid(start: f64, end: f64, n_points: usize) -> Array1<f64> {
      Array1::linspace(start, end, n_points)
    }

    pub(super) fn function(
      range: (f64, f64),
      n_points: usize,
      left: BoundaryBehavior,
      right: BoundaryBehavior,
      neglog: impl Fn(f64) -> f64,
    ) -> Result<Distribution<NegLog>, Report> {
      let y = grid(range.0, range.1, n_points).mapv(neglog);
      let function = DistributionFunction::from_range_values(range, y)?
        .with_left_extrap(left)?
        .with_right_extrap(right)?;
      Ok(Distribution::Function(function))
    }

    pub(super) fn grid_region(distribution: &Distribution<NegLog>, fraction: f64) -> Result<PosteriorRegion, Report> {
      match max_posterior_region(distribution, fraction)? {
        Some(MaxPosteriorRegion::Grid(region)) => Ok(region),
        other => make_error!("Expected a grid region, got {other:?}"),
      }
    }

    pub(super) fn region_errors(
      cells: &[usize],
      distribution: impl Fn(usize) -> Result<Distribution<NegLog>, Report>,
      expected: (f64, f64),
    ) -> Result<Vec<f64>, Report> {
      cells
        .iter()
        .map(|&cells| {
          let region = grid_region(&distribution(cells)?, 0.9)?;
          Ok((region.lower - expected.0).abs().max((region.upper - expected.1).abs()))
        })
        .try_collect()
    }

    pub(super) fn assert_second_order(errors: &[f64]) {
      let coarsest = errors[0].abs();
      let bounded = errors
        .iter()
        .zip(0_i32..)
        .all(|(error, refinements)| error.abs() <= 1.25 * coarsest * 4.0_f64.powi(-refinements));
      let finest_order = (errors[errors.len() - 2] / errors[errors.len() - 1]).abs().log2();
      assert!(
        bounded && finest_order >= 1.9,
        "errors {errors:?} do not fall at second order (finest order {finest_order})"
      );
    }

    pub(super) fn trapezoid(width: f64, y0: f64, y1: f64) -> f64 {
      width * f64::midpoint((-y0).exp(), (-y1).exp())
    }

    pub(super) fn model_density(t: &Array1<f64>, y: &Array1<f64>, x: f64) -> f64 {
      let i = t.iter().position(|&ti| ti > x).unwrap_or(t.len() - 1).max(1);
      let (p0, p1) = ((-y[i - 1]).exp(), (-y[i]).exp());
      p0 + (x - t[i - 1]) / (t[i] - t[i - 1]) * (p1 - p0)
    }
  }
}
