#[cfg(test)]
mod tests {
  use crate::optimize::likelihood::OptimizationMetrics;
  use crate::pretty_assert_ulps_eq;
  use treetime_primitives::LogLh;

  #[test]
  fn test_optimization_metrics_add_sums_each_component_across_partitions() {
    let mut total = OptimizationMetrics::default();

    total.add(&OptimizationMetrics::new(LogLh::new(-10.0), 0.5, -1.0));
    total.add(&OptimizationMetrics::new(LogLh::new(-15.0), 0.3, -2.0));
    total.add(&OptimizationMetrics::new(LogLh::new(-5.0), -0.1, -0.5));

    pretty_assert_ulps_eq!(-10.0 + -15.0 + -5.0, total.log_lh.value(), max_ulps = 4);
    pretty_assert_ulps_eq!(0.5 + 0.3 + -0.1, total.derivative, max_ulps = 4);
    pretty_assert_ulps_eq!(-1.0 + -2.0 + -0.5, total.second_derivative, max_ulps = 4);
  }
}
