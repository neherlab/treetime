#[cfg(test)]
mod tests {
  use crate::distribution::distribution::Distribution;
  use crate::distribution::scaled_distribution::ScaledDistribution;
  use crate::distribution::scaled_distribution_convolution::scaled_distribution_convolution;
  use crate::distribution::y_axis_policy::Plain;
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
  fn test_scaled_distribution_convolve_points() {
    let a = make_point(2.0, 3.0);
    let b = make_point(5.0, 4.0);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert_ulps_eq!(result.peak_value(), 12.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_convolve_functions() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 1.0, 1.0]);
    let b = make_function(array![0.0, 1.0], array![1.0, 1.0]);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }

  #[test]
  fn test_scaled_distribution_convolve_empty() {
    let a = ScaledDistribution::default();
    let b = make_point(1.0, 2.0);
    let result = scaled_distribution_convolution(&a, &b).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_convolve_preserves_normalization() {
    let a = make_function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 10.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![5.0, 20.0, 5.0]);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert_ulps_eq!(result.inner().max_value(), 1.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_convolve_symmetric() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![2.0, 1.0, 2.0]);

    let result_ab = scaled_distribution_convolution(&a, &b).unwrap();
    let result_ba = scaled_distribution_convolution(&b, &a).unwrap();

    assert_ulps_eq!(result_ab.log_scale(), result_ba.log_scale(), max_ulps = 4);
  }
}
