#[cfg(test)]
mod tests {
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{Distribution, DistributionFunction, trim_to_mass};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::{Array1, array, s};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use statrs::distribution::{ContinuousCDF, Normal};
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side, SoftTailLaw};
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_array_eq, pretty_assert_ulps_eq};

  const TRIM_EPS: f64 = 5e-4;

  #[test]
  fn test_mass_trim_cuts_hard_sides_at_the_grid_points_outside_the_quantiles() -> Result<(), Report> {
    let input = helpers::gaussian(-50.0, 10.0, 6001)?;
    let profile = MassProfile::new(&input)?;
    let fractions = profile
      .fractions(TailMass::Required)
      .expect("the input has a finite mass");
    let t = input.t();
    let first = t
      .iter()
      .rposition(|&x| x <= fractions.quantile(TRIM_EPS))
      .expect("a point below the cut");
    let last = t
      .iter()
      .position(|&x| x >= fractions.survival_inverse(TRIM_EPS))
      .expect("a point above the cut");

    let trimmed = helpers::trim(input.clone())?;

    let expected =
      DistributionFunction::from_start_dx_values(t[first], input.dx(), input.y().slice(s![first..=last]).to_owned())?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)?;
    assert_eq!(expected, trimmed);
    assert!(
      fractions.cdf(trimmed.x_min()) <= TRIM_EPS,
      "at most TRIM_EPS is cut below"
    );
    assert!(
      fractions.survival(trimmed.x_max()) <= TRIM_EPS,
      "at most TRIM_EPS is cut above"
    );
    assert!(
      fractions.cdf(t[first + 1]) > TRIM_EPS,
      "the cut is the last grid point below the quantile"
    );
    assert!(
      fractions.survival(t[last - 1]) > TRIM_EPS,
      "the cut is the first grid point above the quantile"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::spacing_0_1(     601)]
  #[case::spacing_0_05(   1201)]
  #[case::spacing_0_025(  2401)]
  #[case::spacing_0_0125( 4801)]
  #[trace]
  fn test_mass_trim_cut_points_lie_within_one_cell_outside_the_normal_quantiles(
    #[case] n_points: usize,
  ) -> Result<(), Report> {
    let normal = Normal::standard();
    let input = helpers::gaussian(-50.0, 10.0, n_points)?;
    let dx = input.dx();

    let trimmed = helpers::trim(input)?;

    let below = (normal.inverse_cdf(TRIM_EPS) - trimmed.x_min()) / dx;
    let above = (trimmed.x_max() - normal.inverse_cdf(1.0 - TRIM_EPS)) / dx;
    for (side, cells) in [("lower", below), ("upper", above)] {
      assert!((0.0..1.0).contains(&cells), "the {side} cut lies {cells} cells outside the normal quantile");
    }
    Ok(())
  }

  #[test]
  fn test_mass_trim_soft_side_with_its_quantile_beyond_the_grid_keeps_its_edge_and_law() -> Result<(), Report> {
    let input = helpers::with_fitted_right_tail(helpers::gaussian(-50.0, 2.0, 5201)?)?;

    let trimmed = helpers::trim(input.clone())?;

    pretty_assert_abs_diff_eq!(input.x_max(), trimmed.x_max(), epsilon = 1e-12);
    assert_eq!(input.right_extrap(), trimmed.right_extrap());
    assert_eq!(BoundaryBehavior::Hard, trimmed.left_extrap());
    Ok(())
  }

  #[test]
  fn test_mass_trim_soft_side_with_its_quantile_inside_the_grid_holds_the_trimmed_mass_in_its_tail()
  -> Result<(), Report> {
    let input = helpers::with_fitted_right_tail(helpers::gaussian(-50.0, 10.0, 6001)?)?;
    let profile = MassProfile::new(&input)?;
    let fractions = profile
      .fractions(TailMass::Required)
      .expect("the input has a finite mass");
    let t = input.t();
    let last = t
      .iter()
      .position(|&x| x >= fractions.survival_inverse(TRIM_EPS))
      .expect("a point above the cut");

    let trimmed = helpers::trim(input.clone())?;

    pretty_assert_abs_diff_eq!(t[last], trimmed.x_max(), epsilon = 1e-12);
    let law = trimmed.right_extrap().soft_law().expect("a soft right side");
    let trimmed_mass = fractions.mass_between(t[last], f64::INFINITY)?;
    pretty_assert_ulps_eq!((-input.y()[last]).exp() / trimmed_mass, law.slope, max_ulps = 8);
    let kept = fractions.total() - fractions.mass_between(f64::NEG_INFINITY, trimmed.x_min())?;
    pretty_assert_ulps_eq!(kept, helpers::total_mass(&trimmed)?, max_ulps = 8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::soft_right(helpers::plateau(Side::Right), 19.0, SoftTailLaw { slope:  0.5 })]
  #[case::soft_left( helpers::plateau(Side::Left),   1.0, SoftTailLaw { slope: -0.5 })]
  #[trace]
  fn test_mass_trim_soft_side_cut_on_a_plateau_gets_a_decaying_tail_with_the_trimmed_mass(
    #[case] input: DistributionFunction<f64, NegLog>,
    #[case] cut: f64,
    #[case] expected: SoftTailLaw,
  ) -> Result<(), Report> {
    let total = 4.5 + 16.5 * (-20.0_f64).exp();
    pretty_assert_ulps_eq!(total, helpers::total_mass(&input)?, max_ulps = 4);

    let trimmed = helpers::trim_at(input, 1e-9)?;

    let (edge, law) = if expected.slope > 0.0 {
      (trimmed.x_max(), trimmed.right_extrap())
    } else {
      (trimmed.x_min(), trimmed.left_extrap())
    };
    pretty_assert_abs_diff_eq!(cut, edge, epsilon = 1e-12);
    let law = law.soft_law().expect("a soft side at the cut");
    pretty_assert_ulps_eq!(expected.slope, law.slope, max_ulps = 8);
    pretty_assert_ulps_eq!(total, helpers::total_mass(&trimmed)?, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_mass_trim_soft_side_cut_at_a_zero_density_point_moves_out_to_a_positive_density() -> Result<(), Report> {
    let y = Array1::from_iter((0..16).map(|i| match i {
      0..=4 => 0.0,
      5 => f64::INFINITY,
      _ => 10.0,
    }));
    let input = DistributionFunction::from_range_values((0.0, 15.0), y)?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }))?;
    let total = 5.0 + 11.0 * (-10.0_f64).exp();
    pretty_assert_ulps_eq!(total, helpers::total_mass(&input)?, max_ulps = 4);

    let trimmed = helpers::trim_at(input, 1e-3)?;

    pretty_assert_abs_diff_eq!(10.0, trimmed.x_max(), epsilon = 1e-12);
    let law = trimmed.right_extrap().soft_law().expect("a soft right side");
    pretty_assert_ulps_eq!(1.0 / 6.0, law.slope, max_ulps = 8);
    pretty_assert_ulps_eq!(total, helpers::total_mass(&trimmed)?, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_mass_trim_hard_approach_side_with_its_quantile_in_the_tail_law_is_not_cut() -> Result<(), Report> {
    let law = HardApproachLaw { t_hard: 0.0, b: 0.0 };
    let input = helpers::flat_on_one_to_two(law)?;
    let total = helpers::total_mass(&input)?;

    let trimmed = helpers::trim(input)?;

    pretty_assert_abs_diff_eq!(1.0, trimmed.x_min(), epsilon = 1e-15);
    assert_eq!(BoundaryBehavior::HardApproach(law), trimmed.left_extrap());
    pretty_assert_abs_diff_eq!(2.0, total, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(total, helpers::total_mass(&trimmed)?, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_mass_trim_hard_approach_side_with_its_quantile_in_the_first_cell_becomes_hard() -> Result<(), Report> {
    let input = helpers::flat_on_one_to_two(HardApproachLaw { t_hard: 0.9999, b: 0.0 })?;

    let trimmed = helpers::trim(input)?;

    pretty_assert_abs_diff_eq!(1.0, trimmed.x_min(), epsilon = 1e-15);
    assert_eq!(BoundaryBehavior::Hard, trimmed.left_extrap());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::fit_points_below_the_cut(   array![30.0, 30.0, 30.0, 0.0, 300.0, 0.0, f64::INFINITY], 1e-9, BoundaryBehavior::Hard)]
  #[case::stored_law_at_the_grid_edge(array![f64::INFINITY, 1.0, 0.0, 0.0, 0.0],             5e-4, BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 2.0 }))]
  #[trace]
  fn test_mass_trim_keeps_the_law_and_mass_of_every_retained_hard_cell(
    #[case] y: Array1<f64>,
    #[case] eps: f64,
    #[case] left: BoundaryBehavior,
  ) -> Result<(), Report> {
    let x_max = (y.len() - 1).to_f64_exact();
    let input = DistributionFunction::from_range_values((0.0, x_max), y)?
      .with_left_extrap(left)?
      .with_right_extrap(BoundaryBehavior::Hard)?;

    let trimmed = helpers::trim_at(input.clone(), eps)?;

    let t = trimmed.t();
    let cells = || t.iter().copied().tuple_windows::<(f64, f64)>();
    let kept = cells()
      .map(|(a, b)| helpers::neg_log_mass_between(&trimmed, a, b))
      .collect::<Result<Array1<f64>, Report>>()?;
    let original = cells()
      .map(|(a, b)| helpers::neg_log_mass_between(&input, a, b))
      .collect::<Result<Array1<f64>, Report>>()?;
    pretty_assert_array_eq!(original, kept);
    let total = helpers::total_mass(&input)?;
    assert!(
      (helpers::total_mass(&trimmed)? - total).abs() <= 2.0 * eps * total,
      "the trim changes the total mass by at most the trimmed mass"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(                Distribution::empty())]
  #[case::point(                Distribution::point(2.0, 0.0))]
  #[case::range(                Distribution::range((1.0, 3.0), 0.0))]
  #[case::without_finite_peak(  helpers::flat(f64::INFINITY, BoundaryBehavior::Hard))]
  #[case::undeclared_tail(      helpers::flat(0.0, BoundaryBehavior::Error))]
  #[case::flat_soft_tail(       helpers::flat(0.0, BoundaryBehavior::Linear(SoftTailLaw { slope: 0.0 })))]
  #[trace]
  fn test_mass_trim_returns_an_input_without_finite_mass_unchanged(#[case] input: Distribution<NegLog>) {
    let trimmed = trim_to_mass(input.clone(), TRIM_EPS).unwrap();
    assert_eq!(input, trimmed);
  }

  mod helpers {
    use super::*;

    pub(super) fn trim(input: DistributionFunction<f64, NegLog>) -> Result<DistributionFunction<f64, NegLog>, Report> {
      trim_at(input, TRIM_EPS)
    }

    pub(super) fn trim_at(
      input: DistributionFunction<f64, NegLog>,
      eps: f64,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      match trim_to_mass(Distribution::Function(input), eps)? {
        Distribution::Function(trimmed) => Ok(trimmed),
        other => Err(eyre::eyre!("trimming a function must give a function, got {other:?}")),
      }
    }

    pub(super) fn neg_log_mass_between(f: &DistributionFunction<f64, NegLog>, a: f64, b: f64) -> Result<f64, Report> {
      let profile = MassProfile::new(f)?;
      Ok(f.amplitude() + profile.y_peak() + profile.neg_log_mass_between(a, b)?)
    }

    pub(super) fn total_mass(f: &DistributionFunction<f64, NegLog>) -> Result<f64, Report> {
      let profile = MassProfile::new(f)?;
      Ok(profile.fractions(TailMass::Required).expect("a finite mass").total())
    }

    pub(super) fn gaussian(
      x_min: f64,
      x_max: f64,
      n_points: usize,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let t = Array1::linspace(x_min, x_max, n_points);
      DistributionFunction::from_start_dx_values(x_min, t[1] - t[0], t.mapv(|t| 0.5 * t * t))?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)
    }

    pub(super) fn with_fitted_right_tail(
      f: DistributionFunction<f64, NegLog>,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let law = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      f.with_right_extrap(BoundaryBehavior::Linear(law))
    }

    pub(super) fn plateau(soft: Side) -> DistributionFunction<f64, NegLog> {
      let (peak, left, right) = match soft {
        Side::Right => (
          0..=4,
          BoundaryBehavior::Hard,
          BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }),
        ),
        Side::Left => (
          16..=20,
          BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }),
          BoundaryBehavior::Hard,
        ),
      };
      let y = Array1::from_iter((0..=20).map(|i| if peak.contains(&i) { 0.0 } else { 20.0 }));
      DistributionFunction::from_range_values((0.0, 20.0), y)
        .and_then(|f| f.with_left_extrap(left))
        .and_then(|f| f.with_right_extrap(right))
        .unwrap()
    }

    pub(super) fn flat_on_one_to_two(left: HardApproachLaw) -> Result<DistributionFunction<f64, NegLog>, Report> {
      DistributionFunction::from_start_dx_values(1.0, 0.5, array![0.0, 0.0, 0.0])?
        .with_left_extrap(BoundaryBehavior::HardApproach(left))?
        .with_right_extrap(BoundaryBehavior::Hard)
    }

    pub(super) fn flat(y: f64, right: BoundaryBehavior) -> Distribution<NegLog> {
      let f = DistributionFunction::from_start_dx_values(0.0, 1.0, array![y, y, y])
        .and_then(|f| f.with_left_extrap(BoundaryBehavior::Hard))
        .and_then(|f| f.with_right_extrap(right))
        .unwrap();
      Distribution::Function(f)
    }
  }
}
