#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_scaled::scaled::ScaledDistribution;
  use crate::policy::Plain;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_scaled_distribution_from_plain_round_trip() {
    let original = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]).unwrap();

    let scaled = ScaledDistribution::from_plain(&original);
    let recovered = scaled.to_plain();

    if let (Distribution::Function(orig), Distribution::Function(rec)) = (&original, &recovered) {
      for i in 0..3 {
        assert_ulps_eq!(orig.y()[i], rec.y()[i], max_ulps = 4);
      }
    } else {
      panic!("Expected Function distributions");
    }
  }

  #[test]
  fn test_scaled_distribution_log_scale() {
    let dist = Distribution::<Plain>::point(1.0, 100.0);
    let scaled = ScaledDistribution::from_plain(&dist);

    assert_ulps_eq!(scaled.log_scale(), 100.0_f64.ln());
    assert_ulps_eq!(scaled.peak_value(), 100.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_inner_is_normalized() {
    let dist = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]).unwrap();
    let scaled = ScaledDistribution::from_plain(&dist);

    assert_ulps_eq!(scaled.inner().max_value(), 1.0);
  }

  #[test]
  fn test_scaled_distribution_empty() {
    let empty = Distribution::<Plain>::Empty;
    let scaled = ScaledDistribution::from_plain(&empty);

    assert!(scaled.is_empty());
    assert!(scaled.log_scale().is_infinite());
  }

  #[test]
  fn test_scaled_distribution_default() {
    let default = ScaledDistribution::default();
    assert!(default.is_empty());
    assert!(default.log_scale().is_infinite());
  }

  #[test]
  fn test_scaled_distribution_renormalize() {
    let dist = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]).unwrap();
    let mut scaled = ScaledDistribution::from_plain(&dist);

    let original_log_scale = scaled.log_scale();
    scaled.renormalize();

    assert_ulps_eq!(scaled.log_scale(), original_log_scale, max_ulps = 4);
    assert_ulps_eq!(scaled.inner().max_value(), 1.0, max_ulps = 4);
  }

  #[test]
  fn test_scaled_distribution_from_parts() {
    let inner = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![0.25, 1.0, 0.5]).unwrap();
    let log_scale = 2.0_f64.ln();
    let scaled = ScaledDistribution::from_parts(log_scale, inner);

    assert_ulps_eq!(scaled.log_scale(), log_scale);
    assert_ulps_eq!(scaled.peak_value(), 2.0, max_ulps = 4);
  }
}
