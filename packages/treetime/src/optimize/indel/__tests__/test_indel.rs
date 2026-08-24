#[cfg(test)]
mod tests {
  use crate::optimize::indel::*;
  use crate::pretty_assert_neg_inf;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_optimize_indel_poisson_zero_rate() {
    let metrics = poisson_indel_log_lh(3, 0.0, 0.1).expect("valid Poisson parameters");
    pretty_assert_neg_inf!(metrics.log_lh.value());
    assert_abs_diff_eq!(metrics.derivative, 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.second_derivative, 0.0, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_indel_poisson_zero_indels() {
    let mu = 5.0;
    let t = 0.1;
    let metrics = poisson_indel_log_lh(0, mu, t).expect("valid Poisson parameters");
    assert_abs_diff_eq!(metrics.log_lh.value(), -mu * t, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.derivative, -mu, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.second_derivative, 0.0, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_indel_poisson_log_lh_value() {
    // k=2, mu=10, t=0.1 => lambda=1.0
    // log P(2|1.0) = 2*ln(1) - 1 - ln(2!) = 0 - 1 - ln(2) = -1 - 0.6931...
    let metrics = poisson_indel_log_lh(2, 10.0, 0.1).expect("valid Poisson parameters");
    let expected_log_lh = -1.0 - 2.0_f64.ln();
    assert_abs_diff_eq!(metrics.log_lh.value(), expected_log_lh, epsilon = 1e-14);
  }

  #[test]
  fn test_optimize_indel_poisson_derivative() {
    let k = 3;
    let mu = 5.0;
    let t = 0.2;
    let metrics = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    // d/dt = k/t - mu = 3/0.2 - 5 = 15 - 5 = 10
    assert_abs_diff_eq!(metrics.derivative, 10.0, epsilon = 1e-13);
    // d2/dt2 = -k/t^2 = -3/0.04 = -75
    // 0.2 is not exactly representable in float64, so t*t has rounding error
    assert_abs_diff_eq!(metrics.second_derivative, -75.0, epsilon = 1e-12);
  }

  #[test]
  fn test_optimize_indel_poisson_mle_at_optimum() {
    // At the maximum likelihood estimate (MLE) t = k/mu, the derivative should be zero
    let k = 5;
    let mu = 10.0;
    let t_mle = k as f64 / mu; // 0.5
    let metrics = poisson_indel_log_lh(k, mu, t_mle).expect("valid Poisson parameters");
    assert_abs_diff_eq!(metrics.derivative, 0.0, epsilon = 1e-14);
  }

  #[test]
  fn test_optimize_indel_poisson_derivative_positive_near_zero() {
    // For k > 0, derivative should be large positive near t=0
    let metrics = poisson_indel_log_lh(1, 5.0, 1e-6).expect("valid Poisson parameters");
    assert!(metrics.derivative > 1e5);
  }

  #[test]
  fn test_optimize_indel_poisson_second_derivative_negative() {
    // Second derivative is always negative when k > 0 (log-concave)
    let metrics = poisson_indel_log_lh(3, 5.0, 0.5).expect("valid Poisson parameters");
    assert!(metrics.second_derivative < 0.0);
  }

  /// Verify statrs ln_factorial agrees with direct computation for small k.
  #[test]
  fn test_optimize_indel_statrs_ln_factorial() {
    assert_abs_diff_eq!(ln_factorial(0), 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(ln_factorial(1), 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(ln_factorial(2), 2.0_f64.ln(), epsilon = 1e-14);
    assert_abs_diff_eq!(ln_factorial(5), 120.0_f64.ln(), epsilon = 1e-14);
    assert_abs_diff_eq!(ln_factorial(10), 3628800.0_f64.ln(), epsilon = 1e-14);
  }
}
