use crate::distribution::distribution_division::distribution_division;
use crate::distribution::scaled_distribution::ScaledDistribution;
use eyre::Report;
use treetime_utils::make_error;

/// Divide one scaled distribution by another.
///
/// Delegates to distribution_division for the inner (normalized) distributions,
/// then combines log_scales (subtraction for division).
pub fn scaled_distribution_division(
  dividend: &ScaledDistribution,
  divisor: &ScaledDistribution,
) -> Result<ScaledDistribution, Report> {
  if dividend.is_empty() {
    return Ok(ScaledDistribution::default());
  }

  if divisor.is_empty() {
    return make_error!("Cannot divide by empty distribution");
  }

  let quotient_inner = distribution_division(dividend.inner(), divisor.inner())?;

  let max_val = quotient_inner.max_value();
  if max_val <= 0.0 || !max_val.is_finite() {
    return Ok(ScaledDistribution::default());
  }

  let combined_log_scale = dividend.log_scale() - divisor.log_scale() + max_val.ln();
  let normalized_quotient = quotient_inner.normalize();

  Ok(ScaledDistribution::from_parts(combined_log_scale, normalized_quotient))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::distribution::distribution::Distribution;
  use crate::distribution::scaled_distribution_multiplication::scaled_distribution_multiplication;
  use crate::distribution::y_axis_policy::Plain;
  use approx::assert_relative_eq;
  use ndarray::{array, Array1};

  fn make_point(t: f64, amplitude: f64) -> ScaledDistribution {
    let dist = Distribution::<Plain>::point(t, amplitude);
    ScaledDistribution::from_plain(&dist)
  }

  fn make_function(x: Array1<f64>, y: Array1<f64>) -> ScaledDistribution {
    let dist = Distribution::<Plain>::function(x, y).unwrap();
    ScaledDistribution::from_plain(&dist)
  }

  #[test]
  fn test_scaled_distribution_divide_empty_by_any() {
    let empty = ScaledDistribution::default();
    let point = make_point(1.0, 2.0);
    let result = scaled_distribution_division(&empty, &point).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_divide_by_empty_fails() {
    let point = make_point(1.0, 2.0);
    let empty = ScaledDistribution::default();
    let result = scaled_distribution_division(&point, &empty);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("empty"));
  }

  #[test]
  fn test_scaled_distribution_divide_functions() {
    let dividend = make_function(array![0.0, 1.0, 2.0], array![2.5, 10.0, 5.0]);
    let divisor = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);

    let result = scaled_distribution_division(&dividend, &divisor).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
    assert_relative_eq!(result.inner().max_value(), 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_divide_inverse_of_multiply() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![2.0, 2.0, 2.0]);

    let product = scaled_distribution_multiplication(&a, &b).unwrap();
    let recovered = scaled_distribution_division(&product, &b).unwrap();

    assert_relative_eq!(recovered.log_scale(), a.log_scale(), epsilon = 0.1);
  }

  #[test]
  fn test_scaled_distribution_divide_preserves_normalization() {
    let dividend = make_function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]);
    let divisor = make_function(array![0.0, 1.0, 2.0], array![2.0, 4.0, 2.0]);

    let result = scaled_distribution_division(&dividend, &divisor).unwrap();

    assert_relative_eq!(result.inner().max_value(), 1.0, epsilon = 1e-10);
  }
}
