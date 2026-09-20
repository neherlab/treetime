#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_scaled::convolve::scaled_distribution_convolution;
  use crate::distribution_scaled::scaled::ScaledDistribution;
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
  fn test_scaled_distribution_convolve_points() {
    let a = make_point(2.0, 3.0);
    let b = make_point(5.0, 4.0);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    assert_ulps_eq!(result.peak_value(), 12.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_convolve_functions() {
    let a = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let b = make_function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]);
    let result = scaled_distribution_convolution(&a, &b).unwrap();

    if let Distribution::Function(f) = result.inner() {
      assert_ulps_eq!(f.grid().x_min(), 0.0, max_ulps = 4);
      assert_ulps_eq!(f.grid().x_max(), 4.0, max_ulps = 4);
      assert_ulps_eq!(result.inner().max_value(), 1.0, max_ulps = 4);
    } else {
      unreachable!("Expected Function distribution");
    }

    assert!(result.log_scale().is_finite());
    assert!(result.log_scale() > 0.0);
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

    if let (Distribution::Function(f_ab), Distribution::Function(f_ba)) = (result_ab.inner(), result_ba.inner()) {
      assert_ulps_eq!(f_ab.grid().x_min(), f_ba.grid().x_min(), max_ulps = 4);
      assert_ulps_eq!(f_ab.grid().x_max(), f_ba.grid().x_max(), max_ulps = 4);
      assert_eq!(f_ab.y().len(), f_ba.y().len());
      approx::assert_ulps_eq!(f_ab.y().as_slice().unwrap(), f_ba.y().as_slice().unwrap(), max_ulps = 4);
    } else {
      unreachable!("Expected Function distributions");
    }
  }
}
