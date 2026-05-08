#[cfg(test)]
mod tests {
  use crate::clock::clock_model::{ClockModel, ClockModelStats, RegressionStats};
  use crate::commands::timetree::output::confidence::{
    date_uncertainty_due_to_rate, determine_rate_std, quantile_to_zscore,
  };
  use approx::assert_relative_eq;
  use ndarray::array;
  use rstest::rstest;
  use treetime_utils::assert_error;

  // Probit function: z = sqrt(2) * erf_inv(2p - 1)
  // Standard values: p=0.025 -> z=-1.959964, p=0.5 -> z=0, p=0.975 -> z=1.959964

  #[rustfmt::skip]
  #[rstest]
  #[case::lower_2_5pct(0.025,  -1.959964)]
  #[case::lower_5pct(  0.05,   -1.644854)]
  #[case::median(       0.5,    0.0      )]
  #[case::upper_95pct(  0.95,   1.644854 )]
  #[case::upper_97_5pct(0.975,  1.959964 )]
  #[case::boundary_zero(0.0,    0.0      )]
  #[case::boundary_one( 1.0,    0.0      )]
  #[trace]
  fn test_quantile_to_zscore(#[case] p: f64, #[case] expected: f64) {
    let z = quantile_to_zscore(p);
    assert_relative_eq!(z, expected, epsilon = 1e-4);
  }

  // Converts [lower_date, center_date, upper_date] + quantile interval to CI.
  // ci_lower = center + z(p_lo) * |lower - center|
  // ci_upper = center + z(p_hi) * |upper - center|

  #[rustfmt::skip]
  #[rstest]
  #[case::symmetric_1sigma(   [9.0, 10.0, 11.0],  (0.025, 0.975), (8.040036, 11.959964))]
  #[case::symmetric_narrow(   [9.5, 10.0, 10.5],  (0.025, 0.975), (9.020018, 10.979982))]
  #[case::asymmetric(         [8.0, 10.0, 11.0],  (0.025, 0.975), (6.080072, 11.959964))]
  #[case::boundary_quantiles( [9.0, 10.0, 11.0],  (0.0,   1.0),   (10.0,     10.0)     )]
  #[case::equal_dates(        [10.0, 10.0, 10.0], (0.025, 0.975), (10.0,     10.0)     )]
  #[trace]
  fn test_date_uncertainty_due_to_rate(
    #[case] dates: [f64; 3],
    #[case] interval: (f64, f64),
    #[case] (expected_lower, expected_upper): (f64, f64),
  ) {
    let (lower, upper) = date_uncertainty_due_to_rate(dates, interval);
    assert_relative_eq!(lower, expected_lower, epsilon = 1e-4);
    assert_relative_eq!(upper, expected_upper, epsilon = 1e-4);
  }

  #[test]
  fn test_determine_rate_std_explicit_clock_std_dev() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(0.001), false, &clock_model).unwrap();
    assert_relative_eq!(result.unwrap(), 0.001);
  }

  #[test]
  fn test_determine_rate_std_rejects_negative() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(-0.001), false, &clock_model);
    assert_error!(result, "--clock-std-dev must be positive, got -0.001");
  }

  #[test]
  fn test_determine_rate_std_rejects_zero() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(Some(0.0), false, &clock_model);
    assert_error!(result, "--clock-std-dev must be positive, got 0");
  }

  #[test]
  fn test_determine_rate_std_none_without_covariation() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(None, false, &clock_model).unwrap();
    assert!(result.is_none());
  }

  #[test]
  fn test_determine_rate_std_from_covariance_matrix() {
    // cov[0,0] = 1e-6, so rate_std = 1e-3
    let clock_model = ClockModel::for_testing_with_stats(
      0.003,
      0.0,
      ClockModelStats::Estimated(RegressionStats {
        chisq: 0.0,
        r_val: 0.9,
        hessian: array![[1.0, 0.0], [0.0, 1.0]],
        cov: array![[1e-6, 0.0], [0.0, 1.0]],
      }),
    );
    let result = determine_rate_std(None, true, &clock_model).unwrap();
    assert_relative_eq!(result.unwrap(), 1e-3, epsilon = 1e-10);
  }

  #[test]
  fn test_determine_rate_std_none_for_fixed_clock_with_covariation() {
    let clock_model = ClockModel::for_testing(0.003, 0.0);
    let result = determine_rate_std(None, true, &clock_model).unwrap();
    assert!(result.is_none());
  }
}
