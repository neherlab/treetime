#[cfg(test)]
mod tests {
  use crate::optimize::indel::poisson_indel_log_lh;
  use crate::pretty_assert_neg_inf;
  use rstest::rstest;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_relative_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_indel_poisson_without_indel_rate_makes_observed_indels_impossible() {
    let metrics = poisson_indel_log_lh(3, 0.0, 0.1).expect("valid Poisson parameters");
    pretty_assert_neg_inf!(metrics.log_lh.value());
    pretty_assert_abs_diff_eq!(0.0, metrics.derivative, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(0.0, metrics.second_derivative, epsilon = 1e-15);
  }

  #[test]
  fn test_indel_poisson_without_indel_rate_and_without_indels_is_certain() {
    let metrics = poisson_indel_log_lh(0, 0.0, 0.1).expect("valid Poisson parameters");
    pretty_assert_abs_diff_eq!(0.0, metrics.log_lh.value(), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(0.0, metrics.derivative, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(0.0, metrics.second_derivative, epsilon = 1e-15);
  }

  #[test]
  fn test_indel_poisson_without_indels_is_the_log_survival_probability() {
    let mu = 5.0;
    let t = 0.1;
    let metrics = poisson_indel_log_lh(0, mu, t).expect("valid Poisson parameters");
    pretty_assert_abs_diff_eq!(-mu * t, metrics.log_lh.value(), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(-mu, metrics.derivative, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(0.0, metrics.second_derivative, epsilon = 1e-15);
  }

  #[test]
  fn test_indel_poisson_slope_grows_as_one_over_t_near_zero() {
    let metrics = poisson_indel_log_lh(1, 5.0, 1e-6).expect("valid Poisson parameters");
    pretty_assert_ulps_eq!(999_995.0, metrics.derivative, max_ulps = 4);
    pretty_assert_ulps_eq!(-1e12, metrics.second_derivative, max_ulps = 4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unit_rate_two_indels(  (  2,  10.0,   0.1), (   -1.6931471805599452,   9.999999999999998, -199.99999999999997))]
  #[case::unit_rate_three_indels((  3,   5.0,   0.2), (   -2.791759469228055,   10.0,                -74.99999999999999))]
  #[case::rate_above_one(        (  3,   5.0,   0.3), (   -2.075364144903562,    5.0,                -33.333333333333336))]
  #[case::underflowing_pmf(      ( 20, 100.0,  20.0), (-1890.317567269912,     -99.0,                 -0.05))]
  #[case::large_count(           (150,   0.5, 100.0), (  -68.21665503520177,     1.0,                 -0.015))]
  #[trace]
  fn test_indel_poisson_matches_the_closed_form_log_pmf_and_its_derivatives(
    #[case] (k, mu, t): (usize, f64, f64),
    #[case] (log_lh, derivative, second_derivative): (f64, f64, f64),
  ) {
    let metrics = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    pretty_assert_relative_eq!(
      log_lh,
      metrics.log_lh.value(),
      max_relative = 1e-14,
      "k ln(mu t) - mu t - ln k! from bc -l at 45 digits on the exact f64 inputs"
    );
    pretty_assert_relative_eq!(
      derivative,
      metrics.derivative,
      max_relative = 1e-15,
      "k / t - mu from bc -l at 45 digits on the exact f64 inputs"
    );
    pretty_assert_relative_eq!(
      second_derivative,
      metrics.second_derivative,
      max_relative = 1e-15,
      "-k / t^2 from bc -l at 45 digits on the exact f64 inputs"
    );
  }
}
