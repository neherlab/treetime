#[cfg(test)]
mod tests {
  use crate::testing::metrics::config::MetricsConfig;
  use crate::testing::metrics::metrics::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_comprehensive_metrics_creation() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected + 0.1;

    let metrics = ValidationMetrics::new(&x, &actual, &expected, 100.0).unwrap();

    assert_eq!(5, metrics.pointwise.total_points);
    assert_eq!(5, metrics.spatial.total_points);

    let expected_r_squared = 1.0 - 0.05 / 2.8;
    assert_ulps_eq!(
      metrics.aggregate.domain_agreement.quality_metrics.r_squared,
      expected_r_squared,
      max_ulps = 4
    );

    assert_ulps_eq!(metrics.distribution.statistics.abs_error_stats.mean, 0.1, max_ulps = 4);
  }

  #[test]
  fn test_metrics_with_custom_config() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = array![1.0, 2.0, 3.0, 2.0, 1.5];

    let mut config = MetricsConfig::default();
    config.distribution.histogram_bins = 20;
    config.spatial.window_half_width = 2;

    let metrics = ValidationMetrics::new_with_config(&x, &actual, &expected, 50.0, &config).unwrap();

    assert_eq!(20, metrics.distribution.histograms.abs_error_histogram.bin_counts.len());
    assert_ulps_eq!(metrics.aggregate.execution_time_ms, 50.0, max_ulps = 4);
    let expected_sliding_max = array![0.0, 0.0, 0.5, 0.5, 0.5];
    assert_ulps_eq!(expected_sliding_max, metrics.spatial.windowed.sliding_max, max_ulps = 0);
  }
}
