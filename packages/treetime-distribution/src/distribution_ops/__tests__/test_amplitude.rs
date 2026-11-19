#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::DistributionFunction;
  use crate::distribution_ops::convolve::distribution_convolution_fine;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use eyre::Report;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, MaxGridPoints};
  use treetime_ops::FftPlans;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_amplitude_repeak_keeps_values_and_moves_the_peak() -> Result<(), Report> {
    let f = Distribution::function(array![0.0, 1.0, 2.0], array![4.0, 1.0, 3.0])?;
    let repeaked = f.repeak()?;
    let Distribution::Function(g) = &repeaked else {
      panic!("expected a function");
    };
    pretty_assert_abs_diff_eq!(&array![3.0, 0.0, 2.0], g.y(), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(1.0, g.amplitude(), epsilon = 1e-15);
    assert_eq!(f, repeaked);
    Ok(())
  }

  #[test]
  fn test_amplitude_normalize_drops_the_level() -> Result<(), Report> {
    let f = DistributionFunction::from_start_dx_values(0.0, 1.0, array![1.0, 0.0, 2.0])?.with_amplitude(5.0);
    let normalized = Distribution::Function(f).normalize()?;
    pretty_assert_abs_diff_eq!(array![1.0, 0.0, 2.0], normalized.y(), epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::function(Distribution::function(array![0.0, 1.0, 2.0], Array1::from_elem(3, f64::INFINITY)).unwrap())]
  #[case::point(   Distribution::point(1.0, f64::INFINITY))]
  #[trace]
  fn test_amplitude_repeak_of_no_probability_is_empty(#[case] dist: Distribution) -> Result<(), Report> {
    assert_eq!(Distribution::Empty, dist.repeak()?);
    Ok(())
  }

  #[test]
  fn test_amplitude_convolution_adds_log_masses() -> Result<(), Report> {
    let a = Distribution::Function(helpers::gaussian(0.0, 0.5)?.with_amplitude(2.0));
    let b = Distribution::Function(helpers::gaussian(1.0, 0.3)?.with_amplitude(-1.5));
    let conv = distribution_convolution_fine(&a, &b, MaxGridPoints::default(), &FftPlans::new())?
      .with_left_extrap(BoundaryBehavior::Hard)
      .with_right_extrap(BoundaryBehavior::Hard);
    pretty_assert_abs_diff_eq!(a.log_mass()? + b.log_mass()?, conv.log_mass()?, epsilon = 1e-15);
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_amplitude_multiplication_adds_values(#[values(-1.0, 0.0, 0.25, 1.0)] t: f64) -> Result<(), Report> {
    let a = Distribution::Function(helpers::gaussian(0.0, 0.5)?.with_amplitude(2.0));
    let b = Distribution::Function(helpers::gaussian(0.5, 0.4)?.with_amplitude(-1.5));
    let product = distribution_multiplication(&a, &b, MaxGridPoints::default())?;
    pretty_assert_abs_diff_eq!(a.eval(t)? + b.eval(t)?, product.eval(t)?, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::point_is_minus_the_amplitude(  Distribution::point(3.0, 0.7),        -0.7)]
  #[case::range_adds_the_log_width(      Distribution::range((1.0, 3.0), 0.7), 2.0_f64.ln() - 0.7)]
  #[trace]
  fn test_amplitude_log_mass_of_point_and_range(
    #[case] dist: Distribution,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    pretty_assert_abs_diff_eq!(expected, dist.log_mass()?, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_amplitude_log_mass_of_gaussian() -> Result<(), Report> {
    let f = Distribution::Function(helpers::gaussian(0.0, 0.5)?.with_amplitude(2.0));
    let expected = (2.0 * std::f64::consts::PI).sqrt().ln() + 0.5_f64.ln() - 2.0;
    pretty_assert_abs_diff_eq!(expected, f.log_mass()?, epsilon = 1e-9);
    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

    const GAUSSIAN_DX: f64 = 0.01;

    pub(super) fn gaussian(mean: f64, sd: f64) -> Result<DistributionFunction, Report> {
      let n = f64_to_usize_saturating((20.0 * sd / GAUSSIAN_DX).round()) + 1;
      let start = mean - 10.0 * sd;
      let y = Array1::from_shape_fn(n, |i| {
        (start + GAUSSIAN_DX * i.to_f64_exact() - mean).powi(2) / (2.0 * sd * sd)
      });
      Ok(
        DistributionFunction::from_start_dx_values(start, GAUSSIAN_DX, y)?
          .with_left_extrap(BoundaryBehavior::Hard)
          .with_right_extrap(BoundaryBehavior::Hard),
      )
    }
  }
}
