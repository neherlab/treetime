#[cfg(test)]
mod tests {
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{BoundaryBehavior, DistributionFunction};
  use approx::assert_relative_eq;
  use eyre::Report;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{HardApproachLaw, Side, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  const HARD: BoundaryBehavior = BoundaryBehavior::Hard;
  const HARD_APPROACH_AT_0: BoundaryBehavior = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 2.0 });

  #[test]
  fn test_mass_profile_neg_log_mass_beyond_error_side_is_error() -> Result<(), Report> {
    let f = helpers::function((0.0, 1.0), array![0.0, 2.0], HARD, BoundaryBehavior::Error)?;
    assert_error!(
      MassProfile::new(&f)?.neg_log_mass_between(0.5, 1.5),
      "Mass between 0.5 and 1.5 reaches above the grid edge 1, where the right side has no tail law"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left_tail_cdf(          (1.0,  2.0 ), HARD_APPROACH_AT_0, HARD,               1.0 )]
  #[case::right_tail_survival(    (-2.0, -1.0), HARD,               HARD_APPROACH_AT_0, -1.0)]
  #[trace]
  fn test_mass_profile_hard_approach_small_targets_are_placed_from_hard_point(
    #[case] range: (f64, f64),
    #[case] left: BoundaryBehavior,
    #[case] right: BoundaryBehavior,
    #[case] outward: f64,
  ) -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::function(range, array![0.0, 0.0], left, right)?)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let total = 1.0 + 1.0 / 3.0;
    let fraction_within = |delta: f64| delta.powi(3) / 3.0 / total;
    let between_1e_17_and_2e_17 = fraction_within(2e-17) - fraction_within(1e-17);
    let (at, from, fraction) = if outward > 0.0 {
      (
        fractions.quantile(fraction_within(1e-17)),
        fractions.quantile_from(1e-17, between_1e_17_and_2e_17),
        fractions.cdf(1e-17),
      )
    } else {
      (
        fractions.survival_inverse(fraction_within(1e-17)),
        fractions.survival_inverse_from(-1e-17, between_1e_17_and_2e_17),
        fractions.survival(-1e-17),
      )
    };
    assert_relative_eq!(outward * 1e-17, at, epsilon = 0.0, max_relative = 1e-14);
    assert_relative_eq!(outward * 2e-17, from, epsilon = 0.0, max_relative = 1e-14);
    assert_relative_eq!(fraction_within(1e-17), fraction, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left_one_ulp_at_edge(        Side::Left,  (1e12 - 0.0001220703125, 1e12                  ))]
  #[case::left_one_ulp_inside_tail(    Side::Left,  (5e11 - 0.00006103515625, 5e11                 ))]
  #[case::left_ten_ulps_at_edge(       Side::Left,  (1e12 - 0.001220703125,  1e12                  ))]
  #[case::right_one_ulp_at_edge(       Side::Right, (-1e12,                  -1e12 + 0.0001220703125))]
  #[case::right_one_ulp_inside_tail(   Side::Right, (-5e11,                  -5e11 + 0.00006103515625))]
  #[case::right_ten_ulps_at_edge(      Side::Right, (-1e12,                  -1e12 + 0.001220703125 ))]
  #[trace]
  fn test_mass_profile_short_interval_in_long_hard_approach_tail_matches_polynomial_integral(
    #[case] side: Side,
    #[case] (a, b): (f64, f64),
  ) -> Result<(), Report> {
    let f = match side {
      Side::Left => helpers::function((1e12, 1e12 + 1.0), array![0.0, 0.0], HARD_APPROACH_AT_0, HARD)?,
      Side::Right => helpers::function((-1e12 - 1.0, -1e12), array![0.0, 0.0], HARD, HARD_APPROACH_AT_0)?,
    };
    let (delta_near, delta_far) = match side {
      Side::Left => (b, a),
      Side::Right => (-a, -b),
    };
    let delta_edge: f64 = 1e12;
    let expected = (delta_near - delta_far)
      * (delta_near.powi(2) + delta_near * delta_far + delta_far.powi(2))
      / (3.0 * delta_edge * delta_edge);
    let profile = MassProfile::new(&f)?;
    assert_relative_eq!(expected, profile.mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    assert_relative_eq!(-expected.ln(), profile.neg_log_mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::whole_support(              1e-200, (0.0,    f64::INFINITY))]
  #[case::grid_only(                  1e-200, (0.0,    1e-150       ))]
  #[case::tail_only(                  1e-200, (1e-150, f64::INFINITY))]
  #[case::finite_tail_part(           1e-200, (1e-150, 1e200        ))]
  #[case::subnormal_slope_whole(      1e-310, (0.0,    f64::INFINITY))]
  #[case::subnormal_slope_tail_part(  1e-310, (1e-150, 1e200        ))]
  #[trace]
  fn test_mass_profile_dominant_tail_beyond_steep_cell_keeps_its_linear_mass(
    #[case] slope: f64,
    #[case] (a, b): (f64, f64),
  ) -> Result<(), Report> {
    let right = BoundaryBehavior::Linear(SoftTailLaw { slope });
    let f = helpers::function((0.0, 1e-150), array![0.0, 800.0], HARD, right)?;
    let profile = MassProfile::new(&f)?;
    let (grid_lo, grid_hi) = (a.max(0.0), b.min(1e-150));
    let density = |t: f64| 1.0 - t / 1e-150;
    let grid = (grid_hi - grid_lo) * f64::midpoint(density(grid_lo), density(grid_hi));
    let (tail_lo, tail_hi) = (a.max(1e-150) - 1e-150, b.max(1e-150) - 1e-150);
    let tail = (-800.0 - slope.ln() - slope * tail_lo).exp() * -(-slope * (tail_hi - tail_lo)).exp_m1();
    let expected = grid + tail;
    let mass = profile.mass_between(a, b)?;
    assert_relative_eq!(expected, mass, epsilon = 0.0, max_relative = 1e-13);
    assert_relative_eq!(mass, (-profile.neg_log_mass_between(a, b)?).exp(), epsilon = 0.0, max_relative = 1e-13);
    if (a, b) == (0.0, f64::INFINITY) {
      let total = profile.fractions(TailMass::Required).unwrap().total();
      assert_relative_eq!(expected, total, epsilon = 0.0, max_relative = 1e-13);
    }
    Ok(())
  }

  #[test]
  fn test_mass_profile_peak_ordinate_1000_gives_relative_masses_of_normalized_copy() -> Result<(), Report> {
    let tails = (
      BoundaryBehavior::Linear(SoftTailLaw { slope: -0.5 }),
      BoundaryBehavior::Linear(SoftTailLaw { slope: 2.0 }),
    );
    let normalized = helpers::function((0.0, 2.0), array![0.5, 0.0, 1.5], tails.0, tails.1)?;
    let shifted = helpers::function((0.0, 2.0), array![1000.5, 1000.0, 1001.5], tails.0, tails.1)?;
    let normalized = MassProfile::new(&normalized)?;
    let shifted = MassProfile::new(&shifted)?;
    let total = normalized.fractions(TailMass::Required).unwrap().total();
    let total_shifted = shifted.fractions(TailMass::Required).unwrap().total();
    pretty_assert_ulps_eq!(total, total_shifted, max_ulps = 0);
    pretty_assert_ulps_eq!(
      normalized.mass_between(-1.0, 1.3)?,
      shifted.mass_between(-1.0, 1.3)?,
      max_ulps = 0
    );
    pretty_assert_ulps_eq!(1000.0, shifted.y_peak(), max_ulps = 0);
    pretty_assert_abs_diff_eq!(
      1000.0 - total.ln(),
      shifted.neg_log_mass(total_shifted),
      epsilon = 1e-15
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::slope_1e_20_over_unit_distance( 1e-20, (1.0, 2.0), 1.0)]
  #[case::slope_0_over_two_units(         0.0,   (1.0, 3.0), 2.0)]
  #[case::slope_0_from_inside_the_grid(   0.0,   (0.5, 2.0), 1.5)]
  #[trace]
  fn test_mass_profile_flat_soft_tail_has_mass_over_finite_distance(
    #[case] slope: f64,
    #[case] (a, b): (f64, f64),
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let right = BoundaryBehavior::Linear(SoftTailLaw { slope });
    let f = helpers::function((0.0, 1.0), array![0.0, 0.0], HARD, right)?;
    pretty_assert_abs_diff_eq!(expected, MassProfile::new(&f)?.mass_between(a, b)?, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::error_side( BoundaryBehavior::Error                            )]
  #[case::flat_side(  BoundaryBehavior::Linear(SoftTailLaw { slope: 0.0 }))]
  #[trace]
  fn test_mass_profile_undefined_tail_mass_follows_policy(#[case] right: BoundaryBehavior) -> Result<(), Report> {
    let f = helpers::function((0.0, 1.0), array![0.0, 2.0], HARD, right)?;
    let profile = MassProfile::new(&f)?;
    assert!(profile.fractions(TailMass::Required).is_none());
    let fractions = profile.fractions(TailMass::ZeroWhereUndefined).unwrap();
    let grid_mass = f64::midpoint(1.0, (-2.0_f64).exp());
    pretty_assert_abs_diff_eq!(grid_mass, fractions.total(), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(0.0, fractions.survival(1.0), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(1.0, fractions.cdf(5.0), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(1.0, fractions.survival_inverse(0.0), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(grid_mass, fractions.mass_between(0.0, 5.0)?, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_mass_profile_interval_beyond_error_side_is_error() -> Result<(), Report> {
    let f = helpers::function((0.0, 1.0), array![0.0, 2.0], BoundaryBehavior::Error, HARD)?;
    let profile = MassProfile::new(&f)?;
    assert_error!(
      profile.mass_between(-0.5, 0.5),
      "Mass between -0.5 and 0.5 reaches below the grid edge 0, where the left side has no tail law"
    );
    let density_at_half = f64::midpoint(1.0, (-2.0_f64).exp());
    pretty_assert_abs_diff_eq!(
      0.5 * f64::midpoint(1.0, density_at_half),
      profile.mass_between(0.0, 0.5)?,
      epsilon = 1e-15
    );
    Ok(())
  }

  #[test]
  fn test_mass_profile_unordered_bounds_are_error() -> Result<(), Report> {
    let f = helpers::function((0.0, 1.0), array![0.0, 2.0], HARD, HARD)?;
    assert_error!(
      MassProfile::new(&f)?.mass_between(0.5, 0.25),
      "Mass between 0.5 and 0.25 requires ordered bounds"
    );
    Ok(())
  }

  #[test]
  fn test_mass_profile_function_without_finite_peak_is_error() -> Result<(), Report> {
    let f = helpers::function((0.0, 1.0), array![f64::INFINITY, f64::INFINITY], HARD, HARD)?;
    assert_error!(
      MassProfile::new(&f),
      "Mass profile requires a finite peak ordinate, got inf"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn function(
      range: (f64, f64),
      y: Array1<f64>,
      left: BoundaryBehavior,
      right: BoundaryBehavior,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      DistributionFunction::from_range_values(range, y)?
        .with_left_extrap(left)?
        .with_right_extrap(right)
    }
  }
}
