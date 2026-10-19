#[cfg(test)]
pub(super) mod tests {
  use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
  use crate::policy::NegLog;
  use crate::{BoundaryBehavior, DistributionFunction};
  use approx::assert_relative_eq;
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::Array1;
  use rstest::rstest;
  use statrs::distribution::{ContinuousCDF, Gamma, Normal};
  use statrs::function::gamma::gamma;
  use std::f64::consts::PI;
  use treetime_grid::{DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};

  const CELLS: [usize; 4] = [100, 200, 400, 800];
  const CELLS_FINE: [usize; 4] = [800, 1600, 3200, 6400];

  #[rustfmt::skip]
  #[rstest]
  #[case::four_hundred_cells(  400 )]
  #[case::eight_hundred_cells( 800 )]
  #[case::fine(                3200)]
  #[trace]
  fn test_mass_profile_gaussian_total_is_exact_to_rounding(#[case] cells: usize) -> Result<(), Report> {
    let f = helpers::sampled((-8.0, 8.0), cells, |t| t * t / 2.0)?;
    let profile = MassProfile::new(&f)?;
    let total = profile.fractions(TailMass::Required).unwrap().total();
    assert_relative_eq!((2.0 * PI).sqrt(), (-profile.neg_log_mass(total)).exp(), epsilon = 0.0, max_relative = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::cdf_below(  -1.0)]
  #[case::cdf_above(  0.5 )]
  #[case::cdf_tail(   2.0 )]
  #[trace]
  fn test_mass_profile_gaussian_cdf_converges_at_second_order(#[case] x: f64) -> Result<(), Report> {
    let normal = Normal::standard();
    let errors = helpers::errors(&CELLS, (-8.0, 8.0), |t| t * t / 2.0, |fractions, _| fractions.cdf(x) - normal.cdf(x))?;
    helpers::assert_second_order(&errors);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unit_rate(  1.0)]
  #[case::fast_rate(  2.5)]
  #[trace]
  fn test_mass_profile_exponential_total_converges_at_second_order(#[case] lambda: f64) -> Result<(), Report> {
    let errors = helpers::errors(&CELLS, (0.0, 30.0), |t| lambda * t - lambda.ln(), |_, total| total - 1.0)?;
    helpers::assert_second_order(&errors);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::total(      None     )]
  #[case::cdf_rising( Some(1.0))]
  #[case::cdf_mode(   Some(2.5))]
  #[case::cdf_tail(   Some(5.0))]
  #[trace]
  fn test_mass_profile_gamma_shape_converges_at_second_order(#[case] x: Option<f64>) -> Result<(), Report> {
    let k = 3.0;
    let gamma_distribution = Gamma::new(k, 1.0)?;
    let errors = helpers::errors(&CELLS_FINE, (0.0, 40.0), |t| t - (k - 1.0) * t.ln(), |fractions, total| match x {
      None => total - gamma(k),
      Some(x) => fractions.cdf(x) - gamma_distribution.cdf(x),
    })?;
    helpers::assert_second_order(&errors);
    Ok(())
  }

  pub(crate) mod helpers {
    use super::*;
    use crate::distribution_ops::mass_profile::MassFractions;

    pub(super) fn errors(
      cells: &[usize],
      range: (f64, f64),
      y: impl Fn(f64) -> f64 + Copy,
      error: impl Fn(&MassFractions<'_>, f64) -> f64,
    ) -> Result<Vec<f64>, Report> {
      cells
        .iter()
        .map(|&cells| {
          let f = sampled(range, cells, y)?;
          let profile = MassProfile::new(&f)?;
          let fractions = profile.fractions(TailMass::Required).unwrap();
          let total = (-profile.neg_log_mass(fractions.total())).exp();
          Ok(error(&fractions, total))
        })
        .try_collect()
    }

    pub(super) fn sampled(
      range: (f64, f64),
      cells: usize,
      y: impl Fn(f64) -> f64,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let t = Array1::linspace(range.0, range.1, cells + 1);
      let f = DistributionFunction::from_range_values(range, t.mapv(y))?;
      let left = if range.0 == 0.0 {
        BoundaryBehavior::Hard
      } else {
        BoundaryBehavior::Linear(SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?)
      };
      let right = BoundaryBehavior::Linear(SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?);
      f.with_left_extrap(left)?.with_right_extrap(right)
    }

    pub(crate) fn assert_second_order(errors: &[f64]) {
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
  }
}
