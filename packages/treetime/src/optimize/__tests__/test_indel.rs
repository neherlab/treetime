#[cfg(test)]
mod tests {
  use crate::optimize::indel::poisson_indel_log_lh;
  use crate::pretty_assert_neg_inf;
  use rstest::rstest;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

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
  #[case::two_indels(  (2, 10.0, 0.1), (-1.0 - 2.0_f64.ln(), 10.0, -200.0))]
  #[case::three_indels((3,  5.0, 0.2), (-1.0 - 6.0_f64.ln(), 10.0,  -75.0))]
  #[trace]
  fn test_indel_poisson_matches_the_closed_form_log_pmf_and_its_derivatives(
    #[case] (k, mu, t): (usize, f64, f64),
    #[case] (log_lh, derivative, second_derivative): (f64, f64, f64),
  ) {
    let metrics = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    pretty_assert_abs_diff_eq!(log_lh, metrics.log_lh.value(), epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(derivative, metrics.derivative, epsilon = 1e-13);
    pretty_assert_abs_diff_eq!(second_derivative, metrics.second_derivative, epsilon = 1e-12);
  }
}
