#[cfg(test)]
mod tests {
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{Distribution, DistributionFunction, trim_to_mass};
  use eyre::Report;
  use ndarray::{Array1, array, s};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use statrs::distribution::{ContinuousCDF, Normal};
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side, SoftTailLaw};
  use treetime_utils::pretty_assert_abs_diff_eq;

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
  fn test_mass_trim_soft_side_with_its_quantile_inside_the_grid_is_cut_and_refit() -> Result<(), Report> {
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

    let trimmed = helpers::trim(input)?;

    pretty_assert_abs_diff_eq!(t[last], trimmed.x_max(), epsilon = 1e-12);
    let refit = SoftTailLaw::fit(trimmed.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
    assert_eq!(BoundaryBehavior::Linear(refit), trimmed.right_extrap());
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
      match trim_to_mass(Distribution::Function(input), TRIM_EPS)? {
        Distribution::Function(trimmed) => Ok(trimmed),
        other => Err(eyre::eyre!("trimming a function must give a function, got {other:?}")),
      }
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
