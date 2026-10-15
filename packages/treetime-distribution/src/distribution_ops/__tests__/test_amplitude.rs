#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::DistributionFunction;
  use crate::distribution_ops::convolve::distribution_convolution_fine;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::policy::NegLog;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array1, array};
  use treetime_grid::{BoundaryBehavior, MaxGridPoints};

  #[test]
  fn test_amplitude_repeak_keeps_values_and_moves_the_peak() -> Result<(), Report> {
    let f = DistributionNegLog::function(array![0.0, 1.0, 2.0], array![4.0, 1.0, 3.0])?;
    let repeaked = f.repeak()?;
    let DistributionNegLog::Function(g) = &repeaked else {
      panic!("expected a function");
    };
    assert_abs_diff_eq!(&array![3.0, 0.0, 2.0], g.y(), epsilon = 1e-15);
    assert_abs_diff_eq!(1.0, g.amplitude(), epsilon = 1e-15);
    assert_eq!(f, repeaked);
    Ok(())
  }

  #[test]
  fn test_amplitude_normalize_drops_the_level() -> Result<(), Report> {
    let f =
      DistributionFunction::<f64, NegLog>::from_start_dx_values(0.0, 1.0, array![1.0, 0.0, 2.0])?.with_amplitude(5.0);
    let normalized = DistributionNegLog::Function(f).normalize()?;
    assert_abs_diff_eq!(array![1.0, 0.0, 2.0], normalized.y()?, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_amplitude_repeak_of_no_probability_is_empty() -> Result<(), Report> {
    let f = DistributionNegLog::function(array![0.0, 1.0, 2.0], Array1::from_elem(3, f64::INFINITY))?;
    assert_eq!(DistributionNegLog::Empty, f.repeak()?);
    assert_eq!(
      DistributionNegLog::Empty,
      DistributionNegLog::point(1.0, f64::INFINITY).repeak()?
    );
    Ok(())
  }

  #[test]
  fn test_amplitude_convolution_adds_log_masses() -> Result<(), Report> {
    let a = DistributionNegLog::Function(helpers::gaussian(0.0, 0.5)?.with_amplitude(2.0));
    let b = DistributionNegLog::Function(helpers::gaussian(1.0, 0.3)?.with_amplitude(-1.5));
    let conv = distribution_convolution_fine(&a, &b, MaxGridPoints::default())?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    // resampling to the finer grid and the FFT tail floor; a lost amplitude would be off by 0.5 or more
    assert_abs_diff_eq!(a.log_mass()? + b.log_mass()?, conv.log_mass()?, epsilon = 1e-4);
    Ok(())
  }

  #[test]
  fn test_amplitude_multiplication_adds_values() -> Result<(), Report> {
    let a = DistributionNegLog::Function(helpers::gaussian(0.0, 0.5)?.with_amplitude(2.0));
    let b = DistributionNegLog::Function(helpers::gaussian(0.5, 0.4)?.with_amplitude(-1.5));
    let product = distribution_multiplication(&a, &b, MaxGridPoints::default())?;
    for t in [-1.0, 0.0, 0.25, 1.0] {
      // linear interpolation on the product grid; a lost amplitude would be off by 1.5 or more
      assert_abs_diff_eq!(a.eval(t)? + b.eval(t)?, product.eval(t)?, epsilon = 1e-3);
    }
    Ok(())
  }

  #[test]
  fn test_amplitude_log_mass_of_point_and_range() -> Result<(), Report> {
    assert_abs_diff_eq!(-0.7, DistributionNegLog::point(3.0, 0.7).log_mass()?, epsilon = 1e-15);
    assert_abs_diff_eq!(
      2.0_f64.ln() - 0.7,
      DistributionNegLog::range((1.0, 3.0), 0.7).log_mass()?,
      epsilon = 1e-15
    );
    Ok(())
  }

  #[test]
  fn test_amplitude_log_mass_of_gaussian() -> Result<(), Report> {
    // peak-1 Gaussian of standard deviation s has mass sqrt(2 pi) s
    let f = DistributionNegLog::Function(helpers::gaussian(0.0, 0.5)?.with_amplitude(2.0));
    let expected = (2.0 * std::f64::consts::PI).sqrt().ln() + 0.5_f64.ln() - 2.0;
    assert_abs_diff_eq!(expected, f.log_mass()?, epsilon = 1e-9);
    Ok(())
  }

  mod helpers {
    use super::*;

    /// Peak-1 Gaussian on mean +- 10 sd (1001 points), hard on both sides
    pub(super) fn gaussian(mean: f64, sd: f64) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let t = Array1::linspace(mean - 10.0 * sd, mean + 10.0 * sd, 1001);
      let y = t.mapv(|t| (t - mean).powi(2) / (2.0 * sd * sd));
      DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)
    }
  }
}
