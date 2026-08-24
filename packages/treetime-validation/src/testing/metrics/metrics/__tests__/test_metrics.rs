#[cfg(test)]
mod tests {
  use crate::testing::metrics::metrics::*;
  use ndarray::array;

  #[test]
  fn test_comprehensive_metrics_creation() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let actual = &expected + 0.1; // Small systematic error

    let metrics = ValidationMetrics::new(&x, &actual, &expected, 100.0).unwrap();

    // Structural checks
    assert_eq!(metrics.pointwise.total_points, 5);
    assert_eq!(metrics.spatial.total_points, 5);

    // R² = 1 - SS_res/SS_tot = 1 - (5*0.1²)/2.8 = 1 - 0.05/2.8
    // mean_expected = 1.8, SS_tot = 0.64+0.04+1.44+0.04+0.64 = 2.8
    let expected_r_squared = 1.0 - 0.05 / 2.8;
    assert_ulps_eq!(
      metrics.aggregate.domain_agreement.quality_metrics.r_squared,
      expected_r_squared,
      max_ulps = 4
    );

    // All absolute errors are 0.1, so mean = 0.1
    assert_ulps_eq!(metrics.distribution.statistics.abs_error_stats.mean, 0.1, max_ulps = 4);
  }

  #[test]
  fn test_metrics_with_custom_config() {
    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 2.0, 1.0];

    let mut config = MetricsConfig::default();
    config.distribution.histogram_bins = 20;
    config.spatial.window_half_width = 2;

    let metrics = ValidationMetrics::new_with_config(&x, &y, &y, 50.0, &config).unwrap();

    assert_eq!(metrics.distribution.histograms.abs_error_histogram.bin_counts.len(), 20);
    assert_ulps_eq!(metrics.aggregate.execution_time_ms, 50.0, max_ulps = 4);
  }
}
