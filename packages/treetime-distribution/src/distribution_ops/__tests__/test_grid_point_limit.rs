#[cfg(test)]
mod tests {
  use crate::distribution_ops::convolve::distribution_convolution_fine;
  use crate::distribution_ops::divide::distribution_division;
  use crate::distribution_ops::mass_domain::resample_to_mass_window;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::distribution_ops::product::distribution_product;
  use eyre::Report;
  use helpers::{flat, flat_function, limit};
  use pretty_assertions::assert_eq;
  use treetime_utils::assert_error;

  const ONE_ABOVE: &str = "A grid over [0, 1000] with spacing 1 needs 1001 points, more than the limit of 1000";

  #[test]
  fn test_grid_point_limit_resample_dx_at_the_limit() -> Result<(), Report> {
    let resampled = flat_function(0.0, 999.0, 2).resample_dx(1.0, limit())?;
    assert_eq!(1_000, resampled.len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_resample_dx_above_the_limit() {
    assert_error!(flat_function(0.0, 1000.0, 2).resample_dx(1.0, limit()), ONE_ABOVE);
  }

  #[test]
  fn test_grid_point_limit_resample_range_dx_at_the_limit() -> Result<(), Report> {
    let resampled = flat_function(0.0, 1000.0, 2).resample_range_dx((0.0, 999.0), 1.0, limit())?;
    assert_eq!(1_000, resampled.len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_resample_range_dx_above_the_limit() {
    assert_error!(
      flat_function(0.0, 1000.0, 2).resample_range_dx((0.0, 1000.0), 1.0, limit()),
      ONE_ABOVE
    );
  }

  #[test]
  fn test_grid_point_limit_convolution_operand_above_the_limit() {
    assert_error!(
      distribution_convolution_fine(&flat(0.0, 999.0, 1_000), &flat(0.0, 1.0, 3), limit()),
      "A grid over [0, 999] with spacing 0.5 needs 1999 points, more than the limit of 1000"
    );
  }

  #[test]
  fn test_grid_point_limit_convolution_result_at_the_limit() -> Result<(), Report> {
    let result = distribution_convolution_fine(&flat(0.0, 499.0, 500), &flat(0.0, 500.0, 501), limit())?;
    assert_eq!(1_000, result.t().len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_convolution_result_above_the_limit() {
    assert_error!(
      distribution_convolution_fine(&flat(0.0, 499.0, 500), &flat(0.0, 501.0, 502), limit()),
      ONE_ABOVE
    );
  }

  #[test]
  fn test_grid_point_limit_multiplication_at_the_limit() -> Result<(), Report> {
    let product = distribution_multiplication(&flat(0.0, 999.0, 1_000), &flat(0.0, 999.0, 1_000), limit())?;
    assert_eq!(1_000, product.t().len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_multiplication_above_the_limit() {
    assert_error!(
      distribution_multiplication(&flat(0.0, 1000.0, 1_001), &flat(0.0, 1000.0, 1_001), limit()),
      ONE_ABOVE
    );
  }

  #[test]
  fn test_grid_point_limit_product_at_the_limit() -> Result<(), Report> {
    let factor = flat(0.0, 999.0, 1_000);
    let product = distribution_product(&[&factor, &factor, &factor], limit())?;
    assert_eq!(1_000, product.t().len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_product_above_the_limit() {
    let factor = flat(0.0, 1000.0, 1_001);
    assert_error!(distribution_product(&[&factor, &factor, &factor], limit()), ONE_ABOVE);
  }

  #[test]
  fn test_grid_point_limit_division_at_the_limit() -> Result<(), Report> {
    let quotient = distribution_division(&flat(0.0, 999.0, 1_000), &flat(0.0, 999.0, 1_000), limit())?;
    assert_eq!(1_000, quotient.t().len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_division_above_the_limit() {
    assert_error!(
      distribution_division(&flat(0.0, 1000.0, 1_001), &flat(0.0, 1000.0, 1_001), limit()),
      ONE_ABOVE
    );
  }

  #[test]
  fn test_grid_point_limit_mass_window_at_the_limit() -> Result<(), Report> {
    let normalized = flat_function(0.0, 1.0, 1_025);
    let windowed = resample_to_mass_window(&normalized, 0.0, 999.0 / 1024.0, 300, limit())?;
    assert_eq!(1_000, windowed.t().len());
    Ok(())
  }

  #[test]
  fn test_grid_point_limit_mass_window_above_the_limit() {
    let normalized = flat_function(0.0, 1.0, 1_025);
    assert_error!(
      resample_to_mass_window(&normalized, 0.0, 1000.0 / 1024.0, 300, limit()),
      "A grid over [0, 0.976563] with spacing 0.000976563 needs 1001 points, more than the limit of 1000"
    );
  }

  mod helpers {
    use crate::__tests__::aliases::DistributionNegLog;
    use crate::distribution_core::function::DistributionFunction;
    use crate::policy::NegLog;
    use ndarray::Array1;
    use treetime_grid::MaxGridPoints;

    pub(super) fn limit() -> MaxGridPoints {
      MaxGridPoints::new(MaxGridPoints::MIN).unwrap()
    }

    pub(super) fn flat(x_min: f64, x_max: f64, n_points: usize) -> DistributionNegLog {
      DistributionNegLog::Function(flat_function(x_min, x_max, n_points))
    }

    pub(super) fn flat_function(x_min: f64, x_max: f64, n_points: usize) -> DistributionFunction<f64, NegLog> {
      DistributionFunction::from_range_values((x_min, x_max), Array1::zeros(n_points)).unwrap()
    }
  }
}
