#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::OptimizationMetrics;
  use crate::pretty_assert_ulps_eq;

  #[test]
  fn test_optimization_metrics_default_is_zero() {
    let metrics = OptimizationMetrics::default();

    pretty_assert_ulps_eq!(metrics.log_lh, 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(metrics.derivative, 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(metrics.second_derivative, 0.0, max_ulps = 4);
  }

  #[test]
  fn test_optimization_metrics_new() {
    let metrics = OptimizationMetrics::new(1.5, -0.3, -2.0);

    pretty_assert_ulps_eq!(metrics.log_lh, 1.5, max_ulps = 4);
    pretty_assert_ulps_eq!(metrics.derivative, -0.3, max_ulps = 4);
    pretty_assert_ulps_eq!(metrics.second_derivative, -2.0, max_ulps = 4);
  }

  #[test]
  fn test_optimization_metrics_add_single() {
    let mut total = OptimizationMetrics::default();
    let other = OptimizationMetrics::new(1.0, 2.0, 3.0);

    total.add(&other);

    pretty_assert_ulps_eq!(total.log_lh, 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(total.derivative, 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(total.second_derivative, 3.0, max_ulps = 4);
  }

  #[test]
  fn test_optimization_metrics_add_multiple_partitions() {
    let mut total = OptimizationMetrics::default();

    let partition1 = OptimizationMetrics::new(-10.0, 0.5, -1.0);
    let partition2 = OptimizationMetrics::new(-15.0, 0.3, -2.0);
    let partition3 = OptimizationMetrics::new(-5.0, -0.1, -0.5);

    total.add(&partition1);
    total.add(&partition2);
    total.add(&partition3);

    let expected_log_lh = -10.0 + -15.0 + -5.0;
    let expected_derivative = 0.5 + 0.3 + -0.1;
    let expected_second_derivative = -1.0 + -2.0 + -0.5;

    pretty_assert_ulps_eq!(total.log_lh, expected_log_lh, max_ulps = 4);
    pretty_assert_ulps_eq!(total.derivative, expected_derivative, max_ulps = 4);
    pretty_assert_ulps_eq!(total.second_derivative, expected_second_derivative, max_ulps = 4);
  }

  #[test]
  fn test_optimization_metrics_add_with_negative_values() {
    let mut total = OptimizationMetrics::new(-100.0, -5.0, -10.0);
    let other = OptimizationMetrics::new(-50.0, 3.0, -5.0);

    total.add(&other);

    pretty_assert_ulps_eq!(total.log_lh, -150.0, max_ulps = 4);
    pretty_assert_ulps_eq!(total.derivative, -2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(total.second_derivative, -15.0, max_ulps = 4);
  }

  #[test]
  fn test_optimization_metrics_add_preserves_precision() {
    let mut total = OptimizationMetrics::default();

    // Add many small values to test precision preservation
    for i in 0..1000 {
      let small = OptimizationMetrics::new(0.001, 0.0001 * (i as f64), -0.00001);
      total.add(&small);
    }

    let expected_log_lh = 1.0;
    let expected_derivative: f64 = (0..1000).map(|i| 0.0001 * (i as f64)).sum();
    let expected_second_derivative = -0.01;

    // Allow more tolerance for 1000 accumulated additions
    pretty_assert_ulps_eq!(total.log_lh, expected_log_lh, max_ulps = 1000);
    pretty_assert_ulps_eq!(total.derivative, expected_derivative, max_ulps = 1000);
    pretty_assert_ulps_eq!(total.second_derivative, expected_second_derivative, max_ulps = 1000);
  }

  #[test]
  fn test_optimization_metrics_add_zero_is_identity() {
    let original = OptimizationMetrics::new(-42.5, 1.23, -4.56);
    let mut total = original.clone();
    let zero = OptimizationMetrics::default();

    total.add(&zero);

    pretty_assert_ulps_eq!(total.log_lh, original.log_lh, max_ulps = 4);
    pretty_assert_ulps_eq!(total.derivative, original.derivative, max_ulps = 4);
    pretty_assert_ulps_eq!(total.second_derivative, original.second_derivative, max_ulps = 4);
  }
}
