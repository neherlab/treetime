#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::scaled_distribution::ScaledDistribution;
  use crate::scaled_distribution_division::scaled_distribution_division;
  use crate::scaled_distribution_multiplication::scaled_distribution_multiplication;
  use crate::y_axis_policy::Plain;
  use approx::assert_ulps_eq;
  use ndarray::{Array1, array};

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
    assert_ulps_eq!(result.inner().max_value(), 1.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_divide_inverse_of_multiply() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![2.0, 2.0, 2.0]);

    let product = scaled_distribution_multiplication(&a, &b).unwrap();
    let recovered = scaled_distribution_division(&product, &b).unwrap();

    assert_ulps_eq!(recovered.log_scale(), a.log_scale(), max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_divide_preserves_normalization() {
    let dividend = make_function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]);
    let divisor = make_function(array![0.0, 1.0, 2.0], array![2.0, 4.0, 2.0]);

    let result = scaled_distribution_division(&dividend, &divisor).unwrap();

    assert_ulps_eq!(result.inner().max_value(), 1.0, max_ulps = 4);
  }
}
