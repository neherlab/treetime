#[cfg(test)]
mod tests {
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{BoundaryBehavior, DistributionFunction};
  use approx::assert_relative_eq;
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{HardApproachLaw, Side, SoftTailLaw};
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  const HARD: BoundaryBehavior = BoundaryBehavior::Hard;

  #[rustfmt::skip]
  #[rstest]
  #[case::rising(         0.7, (2.1,           0.3          ), 0.35 * ((-2.1_f64).exp() + (-0.3_f64).exp()))]
  #[case::falling(        0.7, (0.3,           2.1          ), 0.35 * ((-2.1_f64).exp() + (-0.3_f64).exp()))]
  #[case::flat(           0.5, (1.2,           1.2          ), 0.5 * (-1.2_f64).exp()                       )]
  #[case::steep_rising(   1.0, (1000.0,        0.0          ), 0.5                                          )]
  #[case::steep_falling(  1.0, (0.0,           1000.0       ), 0.5                                          )]
  #[case::infinite_right( 1.0, (0.0,           f64::INFINITY), 1.0                                          )]
  #[case::infinite_left(  1.0, (f64::INFINITY, 0.0          ), 1.0                                          )]
  #[trace]
  fn test_mass_profile_cell_mass_matches_closed_form(
    #[case] h: f64,
    #[case] (y0, y1): (f64, f64),
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let f = helpers::function((0.0, h), array![y0, y1], HARD, HARD)?;
    let profile = MassProfile::new(&f)?;
    let actual = (-profile.neg_log_mass(profile.mass_between(0.0, h)?)).exp();
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::whole_cell(     (0.0,  2.0 ))]
  #[case::left_part(      (0.0,  0.5 ))]
  #[case::right_part(     (1.25, 2.0 ))]
  #[case::inner_part(     (0.3,  1.9 ))]
  #[case::near_left_end(  (1e-9, 2e-9))]
  #[trace]
  fn test_mass_profile_regular_cell_mass_is_integral_of_linear_density(#[case] (a, b): (f64, f64)) -> Result<(), Report> {
    let (y0, y1) = (0.3, 1.7);
    let profile = MassProfile::new(&helpers::function((0.0, 2.0), array![y0, y1], HARD, HARD)?)?;
    let expected = helpers::linear_density_mass(2.0, (y0, y1), (a, b)) / (-y0).exp();
    assert_relative_eq!(expected, profile.mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flat_stored(      0.0, true )]
  #[case::linear_stored(    1.0, true )]
  #[case::fractional_stored(2.5, true )]
  #[case::flat_fitted(      0.0, false)]
  #[case::linear_fitted(    1.0, false)]
  #[case::fractional_fitted(2.5, false)]
  #[trace]
  fn test_mass_profile_hard_cell_holds_mass_of_its_power_law(#[case] beta: f64, #[case] stored: bool) -> Result<(), Report> {
    let anchor_y = 0.4;
    let y = array![anchor_y - beta * 3.0_f64.ln(), anchor_y - beta * 2.0_f64.ln(), anchor_y, f64::INFINITY];
    let right = if stored {
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 3.0, b: beta })
    } else {
      HARD
    };
    let f = helpers::function((0.0, 3.0), y, HARD, right)?;
    let profile = MassProfile::new(&f)?;
    let power = beta + 1.0;
    let below_hard_point = |t: f64| (-anchor_y).exp() * (3.0 - t).powf(power) / power;
    for (a, b) in [(2.0, 3.0), (2.0, 2.5), (2.6, 2.99), (2.999, 3.0)] {
      let expected = (below_hard_point(a) - below_hard_point(b)) / (-profile.y_peak()).exp();
      assert_relative_eq!(expected, profile.mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-14);
    }
    let quadrature = helpers::gauss_legendre(|t| (-(f.interp(t).unwrap() - profile.y_peak())).exp(), (2.0, 3.0), 2000);
    assert_relative_eq!(quadrature, profile.mass_between(2.0, 3.0)?, epsilon = 0.0, max_relative = 1e-13);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flat_near_anchor(        0.0, 1e-10      )]
  #[case::flat_middle(             0.0, 0.5        )]
  #[case::linear_near_hard_point(  1.0, 0.999_999_9)]
  #[case::fractional_middle(       2.5, 0.5        )]
  #[case::fractional_near_anchor(  2.5, 1e-10      )]
  #[case::fractional_near_hard(    2.5, 0.999      )]
  #[trace]
  fn test_mass_profile_hard_cell_fractions_and_inverses_match_power_law(#[case] beta: f64, #[case] x: f64) -> Result<(), Report> {
    let right = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 1.0, b: beta });
    let profile = MassProfile::new(&helpers::function((0.0, 1.0), array![0.0, f64::INFINITY], HARD, right)?)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let power = beta + 1.0;
    let survival = (1.0 - x).powf(power);
    let cdf = -(power * (-x).ln_1p()).exp_m1();
    pretty_assert_abs_diff_eq!(1.0 / power, fractions.total(), epsilon = 1e-15);
    assert_relative_eq!(cdf, fractions.cdf(x), epsilon = 0.0, max_relative = 1e-15);
    assert_relative_eq!(survival, fractions.survival(x), epsilon = 0.0, max_relative = 1e-15);
    if x <= 0.5 {
      assert_relative_eq!(x, fractions.quantile(cdf), epsilon = 0.0, max_relative = 1e-15);
      pretty_assert_abs_diff_eq!(x, fractions.survival_inverse(survival), epsilon = 1e-15);
    } else {
      let distance_to_hard_point = 1.0 - fractions.survival_inverse(survival);
      assert_relative_eq!(1.0 - x, distance_to_hard_point, epsilon = 0.0, max_relative = 1e-15);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::two_points(     2   )]
  #[case::three_points(   3   )]
  #[case::seven_points(   7   )]
  #[case::hundred_cells(  101 )]
  #[case::fine(           4001)]
  #[trace]
  fn test_mass_profile_exponential_total_is_trapezoid_sum_with_exact_tail(#[case] n: usize) -> Result<(), Report> {
    let lambda: f64 = 2.5;
    let t = Array1::linspace(0.0, 3.0, n);
    let y = t.mapv(|ti| -lambda.ln() + lambda * ti);
    let f = DistributionFunction::from_range_values((0.0, 3.0), y)?
      .with_left_extrap(HARD)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: lambda }))?;
    let profile = MassProfile::new(&f)?;
    let total = profile.fractions(TailMass::Required).unwrap().total();
    let h = 3.0 / (n - 1).to_f64_exact();
    let geometric = -(-lambda * h * n.to_f64_exact()).exp_m1() / -(-lambda * h).exp_m1();
    let trapezoid = lambda * h * (geometric - f64::midpoint(1.0, (-3.0 * lambda).exp()));
    let expected = trapezoid + (-3.0 * lambda).exp();
    pretty_assert_abs_diff_eq!(expected, (-profile.neg_log_mass(total)).exp(), epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::whole_support(         (-1.0, f64::INFINITY))]
  #[case::across_all_pieces(     (-0.5, 2.5          ))]
  #[case::inside_power_law(      (-0.8, -0.2         ))]
  #[case::beyond_hard_point(     (-3.0, -1.5         ))]
  #[case::partial_cells(         (0.25, 1.5          ))]
  #[case::inside_one_cell(       (1.25, 1.75         ))]
  #[case::inside_soft_tail(      (2.5,  4.0          ))]
  #[trace]
  fn test_mass_profile_piecewise_model_mass_between_matches_closed_form(#[case] (a, b): (f64, f64)) -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::piecewise_model()?)?;
    pretty_assert_abs_diff_eq!(helpers::piecewise_mass(a, b), profile.mass_between(a, b)?, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::power_law_tail( -0.5, false)]
  #[case::first_cell(     0.25, false)]
  #[case::second_cell(    1.5,  true )]
  #[case::soft_tail(      2.5,  true )]
  #[trace]
  fn test_mass_profile_piecewise_model_fractions_match_closed_form(
    #[case] x: f64,
    #[case] right_of_peak: bool,
  ) -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::piecewise_model()?)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let total = helpers::piecewise_mass(-1.0, f64::INFINITY);
    let cdf = helpers::piecewise_mass(-1.0, x) / total;
    let survival = helpers::piecewise_mass(x, f64::INFINITY) / total;
    pretty_assert_abs_diff_eq!(total, fractions.total(), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(cdf, fractions.cdf(x), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(survival, fractions.survival(x), epsilon = 1e-15);
    let inverse = if right_of_peak { fractions.survival_inverse(survival) } else { fractions.quantile(cdf) };
    pretty_assert_abs_diff_eq!(x, inverse, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::from_power_law_tail( (-0.5, 0.3 ))]
  #[case::from_first_cell(     (0.25, 0.5 ))]
  #[case::into_soft_tail(      (1.5,  0.2 ))]
  #[trace]
  fn test_mass_profile_piecewise_model_quantile_from_lower_bound_holds_fraction(
    #[case] (lower, q): (f64, f64),
  ) -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::piecewise_model()?)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let upper = fractions.quantile_from(lower, q);
    pretty_assert_abs_diff_eq!(q, fractions.mass_between(lower, upper)? / fractions.total(), epsilon = 1e-15);
    let lower_again = fractions.survival_inverse_from(upper, q);
    pretty_assert_abs_diff_eq!(q, fractions.mass_between(lower_again, upper)? / fractions.total(), epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small_cdf_in_rising_cell(       (0.0,  1.0), (1000.0, 0.0   ), false, 0.95  )]
  #[case::small_survival_in_falling_cell( (0.0,  1.0), (0.0,    1000.0), true,  0.05  )]
  #[case::tiny_cdf_in_falling_cell(       (0.0,  1.0), (0.0,    1000.0), false, 1e-23 )]
  #[case::tiny_survival_in_rising_cell(   (-1.0, 0.0), (1000.0, 0.0   ), true,  -1e-23)]
  #[case::tiny_cdf_in_rising_cell(        (0.0,  1.0), (3.0,    0.0   ), false, 1e-23 )]
  #[case::tiny_survival_in_falling_cell(  (-1.0, 0.0), (0.0,    3.0   ), true,  -1e-23)]
  #[trace]
  fn test_mass_profile_small_targets_keep_their_digits(
    #[case] range: (f64, f64),
    #[case] (y0, y1): (f64, f64),
    #[case] from_right: bool,
    #[case] x: f64,
  ) -> Result<(), Report> {
    let f = helpers::function(range, array![y0, y1], HARD, HARD)?;
    let profile = MassProfile::new(&f)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let h = range.1 - range.0;
    let total = helpers::linear_density_mass(h, (y0, y1), (0.0, h));
    let (distance, actual_distance) = if from_right {
      let target = helpers::linear_density_mass_of_right_end(h, (y0, y1), range.1 - x) / total;
      (range.1 - x, range.1 - fractions.survival_inverse(target))
    } else {
      let target = helpers::linear_density_mass(h, (y0, y1), (0.0, x - range.0)) / total;
      (x - range.0, fractions.quantile(target) - range.0)
    };
    assert_relative_eq!(distance, actual_distance, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::falling_near_left_end(  (0.0, 3.0), 1e-10,       false)]
  #[case::falling_middle_cdf(     (0.0, 3.0), 0.5,         false)]
  #[case::falling_middle_survival((0.0, 3.0), 0.5,         true )]
  #[case::falling_near_right_end( (0.0, 3.0), 0.999_999_9, true )]
  #[case::rising_near_left_end(   (3.0, 0.0), 1e-10,       false)]
  #[case::rising_middle_cdf(      (3.0, 0.0), 0.5,         false)]
  #[case::rising_middle_survival( (3.0, 0.0), 0.5,         true )]
  #[case::rising_near_right_end(  (3.0, 0.0), 0.999_999_9, true )]
  #[trace]
  fn test_mass_profile_inverse_is_accurate_near_cell_ends(
    #[case] (y0, y1): (f64, f64),
    #[case] u: f64,
    #[case] from_right: bool,
  ) -> Result<(), Report> {
    let f = helpers::function((0.0, 1.0), array![y0, y1], HARD, HARD)?;
    let profile = MassProfile::new(&f)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let total = helpers::linear_density_mass(1.0, (y0, y1), (0.0, 1.0));
    let inverse = if from_right {
      fractions.survival_inverse(helpers::linear_density_mass(1.0, (y0, y1), (u, 1.0)) / total)
    } else {
      fractions.quantile(helpers::linear_density_mass(1.0, (y0, y1), (0.0, u)) / total)
    };
    pretty_assert_abs_diff_eq!(u, inverse, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rising_left_end(   -1e-12, 1e-10      )]
  #[case::rising_middle(     -1e-12, 0.5        )]
  #[case::falling_left_end(  1e-12,  1e-10      )]
  #[case::falling_middle(    1e-12,  0.5        )]
  #[case::falling_right_end( 1e-12,  0.999_999_9)]
  #[trace]
  fn test_mass_profile_inverse_is_accurate_for_nearly_flat_cells(#[case] slope: f64, #[case] u: f64) -> Result<(), Report> {
    let y0 = if slope < 0.0 { -slope } else { 0.0 };
    let ys = (y0, y0 + slope);
    let f = helpers::function((0.0, 1.0), array![ys.0, ys.1], HARD, HARD)?;
    let profile = MassProfile::new(&f)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let cdf = helpers::linear_density_mass(1.0, ys, (0.0, u)) / helpers::linear_density_mass(1.0, ys, (0.0, 1.0));
    pretty_assert_abs_diff_eq!(u, fractions.quantile(cdf), epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::whole_valley_cells(   (2.0, 5.0), 3.0)]
  #[case::partial_valley_cells( (2.5, 4.5), 2.0)]
  #[case::one_valley_cell(      (3.0, 4.0), 1.0)]
  #[trace]
  fn test_mass_profile_valley_between_peaks_keeps_its_mass(
    #[case] (a, b): (f64, f64),
    #[case] valley_width: f64,
  ) -> Result<(), Report> {
    let y = array![0.0, 0.0, 100.0, 100.0, 100.0, 100.0, 0.0, 0.0];
    let profile = MassProfile::new(&helpers::function((0.0, 7.0), y, HARD, HARD)?)?;
    assert_relative_eq!(valley_width * (-100.0_f64).exp(), profile.mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    assert_relative_eq!(100.0 - valley_width.ln(), profile.neg_log_mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::soft_tail_beyond_underflow(  (1000.0, 1001.0))]
  #[case::cells_beyond_underflow(      (800.0,  1000.0))]
  #[case::cells_and_soft_tail(         (800.0,  1001.0))]
  #[case::edge_cell_and_soft_tail(     (999.0,  1010.0))]
  #[case::near_peak(                   (0.0,    1.0   ))]
  #[case::across_moderate_densities(   (0.0,    700.0 ))]
  #[case::from_inside_first_cell(      (0.5,    700.0 ))]
  #[trace]
  fn test_mass_profile_neg_log_mass_of_sampled_unit_exponential_matches_closed_form(
    #[case] (a, b): (f64, f64),
  ) -> Result<(), Report> {
    let right = BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 });
    let t = Array1::linspace(0.0, 1000.0, 1001);
    let profile = MassProfile::new(&helpers::function((0.0, 1000.0), t, HARD, right)?)?;
    let e = (-1.0_f64).exp();
    let whole_cells = |lo: f64, hi: f64| (lo, f64::midpoint(1.0, e) * -(lo - hi).exp_m1() / (1.0 - e));
    let first = a.floor();
    let mut pieces = vec![];
    if a < 1000.0 && a > first {
      let density_a = 1.0 + (e - 1.0) * (a - first);
      pieces.push((first, (first + 1.0 - a) * f64::midpoint(density_a, e)));
    }
    if a.ceil() < 1000.0 {
      pieces.push(whole_cells(a.ceil(), b.min(1000.0)));
    }
    if b > 1000.0 {
      let start = a.max(1000.0);
      pieces.push((start, -(start - b).exp_m1()));
    }
    let r_min = pieces.iter().map(|(r, _)| *r).fold(f64::INFINITY, f64::min);
    let expected = r_min - pieces.iter().map(|(r, factor)| factor * (r_min - r).exp()).sum::<f64>().ln();
    assert_relative_eq!(expected, profile.neg_log_mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::near_peak(           (0.0,  1.0  ))]
  #[case::across_power_law(    (-0.8, 0.5  ))]
  #[case::across_soft_tail(    (1.5,  4.0  ))]
  #[case::whole_finite_span(   (-1.0, 10.0 ))]
  #[trace]
  fn test_mass_profile_neg_log_mass_agrees_with_linear_mass_near_peak(#[case] (a, b): (f64, f64)) -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::piecewise_model()?)?;
    assert_relative_eq!(-profile.mass_between(a, b)?.ln(), profile.neg_log_mass_between(a, b)?, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[test]
  fn test_mass_profile_neg_log_mass_of_empty_interval_is_infinite() -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::piecewise_model()?)?;
    pretty_assert_ulps_eq!(f64::INFINITY, profile.neg_log_mass_between(-3.0, -1.5)?, max_ulps = 0);
    pretty_assert_ulps_eq!(f64::INFINITY, profile.neg_log_mass_between(0.5, 0.5)?, max_ulps = 0);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::steep_narrow_cell(   999.0, 1e-20,  5e-21,   746.106_637_698_808_4)]
  #[case::wide_cell_tiny_start(  1.0, 1e200,  2e-200, 1_380.857_908_615_867_5)]
  #[case::flat_wide_cell(        0.0, 1e200,  2e-200,  459.823_871_418_249_2)]
  #[trace]
  fn test_mass_profile_hard_cell_mass_far_below_plain_scale_keeps_its_log_mass(
    #[case] beta: f64,
    #[case] h: f64,
    #[case] distance: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let left = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: beta });
    let profile = MassProfile::new(&helpers::function((0.0, h), array![f64::INFINITY, 0.0], left, HARD)?)?;
    assert_relative_eq!(expected, profile.neg_log_mass_between(0.0, distance)?, epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::infinite_first(  (2000.0, 2000.2), array![f64::INFINITY, 3.0, 4.0])]
  #[case::finite(          (2000.0, 2000.2), array![5.0, 3.0, 4.0])]
  #[case::negative_tenths( (-7.3,   -6.6  ), array![f64::INFINITY, 1.0, 2.0, 9.0, f64::INFINITY, 5.0, 6.0, 3.0])]
  #[case::thirds(          (1999.0, 2001.0), array![f64::INFINITY, 1.0, 2.0, 3.0, 4.0, 5.0, f64::INFINITY])]
  #[trace]
  fn test_mass_profile_mass_between_stored_coordinates_sums_whole_cells(
    #[case] range: (f64, f64),
    #[case] y: Array1<f64>,
  ) -> Result<(), Report> {
    let f = helpers::function(range, y.clone(), HARD, HARD)?;
    let profile = MassProfile::new(&f)?;
    let x = f.grid_fn().x();
    let y_peak = y.iter().copied().fold(f64::INFINITY, f64::min);
    let cells: Vec<f64> = (0..y.len() - 1)
      .map(|c| match f.grid_fn().hard_cell(c, HARD, HARD)? {
        Some(hard) => Ok((y_peak - hard.anchor.y).exp() * f.dx() / (hard.law.b + 1.0)),
        None => Ok(helpers::linear_density_mass(f.dx(), (y[c] - y_peak, y[c + 1] - y_peak), (0.0, f.dx()))),
      })
      .try_collect::<_, _, Report>()?;
    for (i, k) in (0..x.len()).tuple_combinations() {
      let expected: f64 = cells[i..k].iter().sum();
      assert_relative_eq!(expected, profile.mass_between(x[i], x[k])?, epsilon = 0.0, max_relative = 1e-15);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beyond_hard_point(        -1.5,  Side::Right, 0.0                                         )]
  #[case::power_law_tail(           -0.5,  Side::Right, (-0.5_f64).exp() * 0.25                     )]
  #[case::left_edge_from_tail(       0.0,  Side::Left,  (-0.5_f64).exp()                            )]
  #[case::left_edge_from_grid(       0.0,  Side::Right, (-0.5_f64).exp()                            )]
  #[case::rising_cell_middle(        0.5,  Side::Left,  f64::midpoint((-0.5_f64).exp(), 1.0)        )]
  #[case::peak(                      1.0,  Side::Left,  1.0                                         )]
  #[case::falling_cell_quarter(      1.25, Side::Right, 1.0 + (-1.5_f64).exp_m1() * 0.25             )]
  #[case::right_edge(                2.0,  Side::Right, (-1.5_f64).exp()                            )]
  #[case::soft_tail(                 3.0,  Side::Left,  (-1.5_f64).exp() * (-0.8_f64).exp()         )]
  #[trace]
  fn test_mass_profile_density_follows_the_cell_and_tail_laws(
    #[case] x: f64,
    #[case] side: Side,
    #[case] density: f64,
  ) -> Result<(), Report> {
    let profile = MassProfile::new(&helpers::piecewise_model()?)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let total = helpers::piecewise_mass(-1.0, f64::INFINITY);
    pretty_assert_ulps_eq!(density / total, fractions.density(x, side), max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::outside_left(  0.0, Side::Left,  0.0)]
  #[case::at_hard_point( 0.0, Side::Right, 1.0)]
  #[case::inside(        0.5, Side::Left,  1.0)]
  #[case::inside_right(  1.0, Side::Left,  1.0)]
  #[case::outside_right( 1.0, Side::Right, 0.0)]
  #[trace]
  fn test_mass_profile_density_of_a_flat_hard_cell_is_one_sided_at_its_ends(
    #[case] x: f64,
    #[case] side: Side,
    #[case] density: f64,
  ) -> Result<(), Report> {
    let flat_at_hard_point = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 0.0 });
    let f = helpers::function((0.0, 1.0), array![f64::INFINITY, 0.0], flat_at_hard_point, HARD)?;
    let profile = MassProfile::new(&f)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    pretty_assert_ulps_eq!(density, fractions.density(x, side), max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flat_earlier_outside(       (0.0, Side::Left ), (0.0, Side::Left ), 0.0                )]
  #[case::flat_earlier_at_hard_point( (0.0, Side::Left ), (0.0, Side::Right), 1.0                )]
  #[case::flat_earlier_inside(        (0.0, Side::Left ), (0.5, Side::Left ), 1.0                )]
  #[case::flat_later_outside(         (0.0, Side::Right), (3.0, Side::Right), 0.0                )]
  #[case::flat_later_at_hard_point(   (0.0, Side::Right), (3.0, Side::Left ), 1.0                )]
  #[case::flat_later_inside(          (0.0, Side::Right), (2.5, Side::Right), 1.0                )]
  #[case::sublinear_at_hard_point(    (0.5, Side::Left ), (0.0, Side::Right), 0.0                )]
  #[case::sublinear_inside(           (0.5, Side::Left ), (0.5, Side::Left ), 0.5_f64.sqrt()     )]
  #[case::quadratic_at_hard_point(    (2.0, Side::Right), (3.0, Side::Left ), 0.0                )]
  #[case::quadratic_inside(           (2.0, Side::Right), (2.5, Side::Right), 0.25               )]
  #[trace]
  fn test_mass_profile_density_of_a_hard_approach_tail_is_its_one_sided_limit_at_the_hard_point(
    #[case] (b, tail_side): (f64, Side),
    #[case] (x, side): (f64, Side),
    #[case] density: f64,
  ) -> Result<(), Report> {
    let (left, right) = match tail_side {
      Side::Left => (BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b }), HARD),
      Side::Right => (HARD, BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 3.0, b })),
    };
    let f = helpers::function((1.0, 2.0), array![0.0, 0.0], left, right)?;
    let profile = MassProfile::new(&f)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();
    let total = 1.0 + 1.0 / (b + 1.0);
    pretty_assert_ulps_eq!(density / total, fractions.density(x, side), max_ulps = 4);
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

    pub(super) fn linear_density_mass(h: f64, (ya, yb): (f64, f64), (lo, hi): (f64, f64)) -> f64 {
      let (pa, pb) = ((-ya).exp(), (-yb).exp());
      let density = |u: f64| {
        if u <= h / 2.0 {
          pa + (pb - pa) * (u / h)
        } else {
          pb + (pa - pb) * ((h - u) / h)
        }
      };
      (hi - lo) * f64::midpoint(density(lo), density(hi))
    }

    pub(super) fn linear_density_mass_of_right_end(h: f64, (ya, yb): (f64, f64), distance: f64) -> f64 {
      let (pa, pb) = ((-ya).exp(), (-yb).exp());
      distance * f64::midpoint(pb + (pa - pb) * (distance / h), pb)
    }

    pub(super) fn gauss_legendre(f: impl Fn(f64) -> f64, (a, b): (f64, f64), panels: usize) -> f64 {
      let nodes = [
        (-(0.6_f64).sqrt(), 5.0 / 9.0),
        (0.0, 8.0 / 9.0),
        ((0.6_f64).sqrt(), 5.0 / 9.0),
      ];
      let width = (b - a) / panels.to_f64_exact();
      (0..panels)
        .map(|i| {
          let center = a + width * (i.to_f64_exact() + 0.5);
          nodes
            .iter()
            .map(|(node, weight)| weight * f(center + node * width / 2.0))
            .sum::<f64>()
            * width
            / 2.0
        })
        .sum()
    }

    pub(super) fn piecewise_model() -> Result<DistributionFunction<f64, NegLog>, Report> {
      function(
        (0.0, 2.0),
        array![0.5, 0.0, 1.5],
        BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: -1.0, b: 2.0 }),
        BoundaryBehavior::Linear(SoftTailLaw { slope: 0.8 }),
      )
    }

    pub(super) fn piecewise_mass(a: f64, b: f64) -> f64 {
      let power_law = |t: f64| (-0.5_f64).exp() * (t + 1.0).clamp(0.0, 1.0).powi(3) / 3.0;
      let left = power_law(b.min(0.0)) - power_law(a.min(0.0));
      let cell = |offset: f64, ys: (f64, f64), lo: f64, hi: f64| {
        let (lo, hi) = (lo.max(offset), hi.min(offset + 1.0));
        if lo < hi {
          linear_density_mass(1.0, ys, (lo - offset, hi - offset))
        } else {
          0.0
        }
      };
      let grid = cell(0.0, (0.5, 0.0), a, b) + cell(1.0, (0.0, 1.5), a, b);
      let soft = |t: f64| (-1.5_f64).exp() * (-0.8 * (t - 2.0).max(0.0)).exp() / 0.8;
      let right = soft(a.max(2.0)) - soft(b.max(2.0));
      left + grid + right
    }
  }
}
