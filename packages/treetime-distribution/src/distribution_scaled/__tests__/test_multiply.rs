#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_scaled::distribution_scaled::ScaledDistribution;
  use crate::distribution_scaled::multiply::{scaled_distribution_multiplication, scaled_distribution_multiply_many};
  use crate::policy::Plain;
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
  fn test_scaled_distribution_multiply_points_same_location() {
    let a = make_point(2.0, 3.0);
    let b = make_point(2.0, 4.0);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert_ulps_eq!(result.peak_value(), 12.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_multiply_points_different_location() {
    let a = make_point(2.0, 3.0);
    let b = make_point(5.0, 4.0);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_multiply_functions() {
    // [1,2,1] * [2,1,2] = [2,2,2] (elementwise), normalized to [1,1,1]
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![2.0, 1.0, 2.0]);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    // Verify actual y-values after normalization
    // Normalized inputs: a_norm = [0.5, 1, 0.5], b_norm = [1, 0.5, 1]
    // Product: [0.5, 0.5, 0.5], max=0.5, normalized: [1, 1, 1]
    if let Distribution::Function(f) = result.inner() {
      assert_ulps_eq!(f.y()[0], 1.0, max_ulps = 4);
      assert_ulps_eq!(f.y()[1], 1.0, max_ulps = 4);
      assert_ulps_eq!(f.y()[2], 1.0, max_ulps = 4);
    } else {
      unreachable!("Expected Function distribution");
    }

    // Verify log_scale is finite and produces correct peak value
    // Peak of product = 2 * 2 * 0.5 = 2 (from a.max * b.max * normalized_product_max)
    assert!(result.log_scale().is_finite());
    assert_ulps_eq!(result.peak_value(), 2.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_multiply_many_no_underflow() {
    let dists: Vec<ScaledDistribution> = std::iter::repeat_with(|| make_point(0.0, 0.0001)).take(100).collect();
    let refs: Vec<&ScaledDistribution> = dists.iter().collect();

    let result = scaled_distribution_multiply_many(&refs).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());

    let expected_log_scale = 100.0 * 0.0001_f64.ln();
    assert_ulps_eq!(result.log_scale(), expected_log_scale, max_ulps = 16);
  }

  #[test]
  fn test_scaled_distribution_multiply_empty_propagates() {
    let a = ScaledDistribution::default();
    let b = make_point(1.0, 2.0);
    let result = scaled_distribution_multiplication(&a, &b).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_multiply_preserves_normalization() {
    let a = make_function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![5.0, 10.0, 5.0]);

    let result = scaled_distribution_multiplication(&a, &b).unwrap();

    assert_ulps_eq!(result.inner().max_value(), 1.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_multiply_many_empty() {
    let result = scaled_distribution_multiply_many(&[]).unwrap();
    assert!(result.is_empty());
  }

  #[test]
  fn test_scaled_distribution_multiply_many_single() {
    let a = make_point(1.0, 5.0);
    let result = scaled_distribution_multiply_many(&[&a]).unwrap();

    assert_ulps_eq!(result.log_scale(), a.log_scale(), max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_multiply_many_functions_same_grid() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![0.5, 1.0, 0.5]);
    let c = make_function(array![0.0, 1.0, 2.0], array![0.5, 1.0, 0.5]);

    let result = scaled_distribution_multiply_many(&[&a, &b, &c]).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
    assert_ulps_eq!(result.inner().max_value(), 1.0, max_ulps = 4);

    if let Distribution::Function(f) = result.inner() {
      assert_ulps_eq!(f.y()[1], 1.0, max_ulps = 4);
      assert_ulps_eq!(f.y()[0], 0.125, max_ulps = 4);
      assert_ulps_eq!(f.y()[2], 0.125, max_ulps = 4);
    } else {
      unreachable!("Expected Function distribution");
    }
  }

  #[test]
  fn test_scaled_distribution_multiply_many_functions_underflow_resistance() {
    let small_peak = 1e-50;
    let dists: Vec<ScaledDistribution> = std::iter::repeat_with(|| {
      make_function(
        array![0.0, 1.0, 2.0],
        array![small_peak * 0.5, small_peak, small_peak * 0.5],
      )
    })
    .take(10)
    .collect();
    let refs: Vec<&ScaledDistribution> = dists.iter().collect();

    let result = scaled_distribution_multiply_many(&refs).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());

    let expected_log_scale = 10.0 * small_peak.ln();
    assert_ulps_eq!(result.log_scale(), expected_log_scale, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_multiply_many_mixed_types_fallback() {
    let func = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let point = make_point(1.0, 3.0);

    let result = scaled_distribution_multiply_many(&[&func, &point]).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }

  #[test]
  fn test_scaled_distribution_multiply_many_different_grids_fallback() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 0.5, 1.0, 1.5, 2.0], array![1.0, 1.5, 2.0, 1.5, 1.0]);

    let result = scaled_distribution_multiply_many(&[&a, &b]).unwrap();

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
  }
}
