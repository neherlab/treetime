#[cfg(test)]
mod tests {

  use crate::optimize::indel::poisson_indel_log_lh;
  use crate::optimize::likelihood::evaluate_partitions;
  use crate::optimize::method_newton::chain_rule_sqrt;
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use ndarray::array;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::{
    pretty_assert_abs_diff_eq, pretty_assert_relative_eq, pretty_assert_ulps_eq, prop_assert_relative_eq,
    prop_assert_ulps_eq,
  };

  #[test]
  fn test_optimize_method_chain_rule_at_zero() {
    let (ds, d2s) = chain_rule_sqrt(0.0, 100.0, -500.0);
    pretty_assert_abs_diff_eq!(0.0, ds, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(200.0, d2s, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_method_chain_rule_analytical() {
    let (ds, d2s) = chain_rule_sqrt(0.3, 10.0, -100.0);
    pretty_assert_abs_diff_eq!(6.0, ds, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(-16.0, d2s, epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_matches_closed_form_in_sqrt_t(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::mixture(
      coefficients,
      helpers::jc69_rate_spectrum(),
    ));
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_partitions(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let expected = helpers::closed_form_in_sqrt_t(s, (0.5, 0.5), k, mu);
    pretty_assert_relative_eq!(
      expected.first,
      dl_ds,
      max_relative = 1e-14,
      "dL/ds against the derivative in s of ln(a + b exp(-4/3 s^2)) + ln Poisson(k; mu s^2)"
    );
    pretty_assert_relative_eq!(
      expected.second,
      d2l_ds2,
      max_relative = 1e-14,
      "d2L/ds2 against the second derivative in s of ln(a + b exp(-4/3 s^2)) + ln Poisson(k; mu s^2)"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::second_cancels_to_zero(    (0.5, -3.0,                     4.0), (-2.0,                     0.0              ))]
  #[case::second_cancels_to_residual((0.5, -3.0 + 2.0_f64.powi(-40), 4.0), (-2.0 + 2.0_f64.powi(-40), 2.0_f64.powi(-39)))]
  #[case::first_cancels_to_zero(     (0.5, -1.0,                     4.0), ( 0.0,                     4.0              ))]
  #[trace]
  fn test_optimize_method_chain_rule_sqrt_exact_when_mixed_signs_cancel(
    #[case] (s, alpha, gamma): (f64, f64, f64),
    #[case] (expected_dl_ds, expected_d2l_ds2): (f64, f64),
  ) {
    let t = s * s;
    let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, alpha + gamma * t, gamma);
    pretty_assert_ulps_eq!(
      expected_dl_ds,
      dl_ds,
      max_ulps = 0,
      "2 alpha s + 2 gamma s^3, the direct s-derivative of alpha s^2 + gamma s^4 / 2, exact for dyadic inputs"
    );
    pretty_assert_ulps_eq!(
      expected_d2l_ds2,
      d2l_ds2,
      max_ulps = 0,
      "2 alpha + 6 gamma s^2, the direct second s-derivative of alpha s^2 + gamma s^4 / 2, exact for dyadic inputs"
    );
  }

  proptest! {
    #[test]
    fn test_prop_optimize_method_chain_rule_sqrt_matches_direct_derivative_of_quadratic(
      s in generators::gen_s(),
      (alpha, gamma) in generators::gen_same_sign_quadratic(),
    ) {
      let t = s * s;
      let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, alpha + gamma * t, gamma);
      let direct_dl_ds = 2.0 * alpha * s + 2.0 * gamma * s * s * s;
      let direct_d2l_ds2 = 2.0 * alpha + 6.0 * gamma * s * s;
      prop_assert_relative_eq!(
        direct_dl_ds,
        dl_ds,
        max_relative = 1e-15,
        "2 alpha s + 2 gamma s^3, the direct s-derivative of alpha s^2 + gamma s^4 / 2"
      );
      prop_assert_relative_eq!(
        direct_d2l_ds2,
        d2l_ds2,
        max_relative = 1e-15,
        "2 alpha + 6 gamma s^2, the direct second s-derivative of alpha s^2 + gamma s^4 / 2"
      );
    }

    #[test]
    fn test_prop_optimize_method_chain_rule_sqrt_linear_under_power_of_two_scaling(
      s in generators::gen_s(),
      dl_dt in generators::gen_dl_dt(),
      d2l_dt2 in generators::gen_d2l_dt2(),
      k in generators::gen_power_of_two(),
    ) {
      let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
      let (dl_ds_k, d2l_ds2_k) = chain_rule_sqrt(s, k * dl_dt, k * d2l_dt2);
      prop_assert_ulps_eq!(
        k * dl_ds,
        dl_ds_k,
        max_ulps = 0,
        "scaling by a power of two is exact while every intermediate stays normal"
      );
      prop_assert_ulps_eq!(
        k * d2l_ds2,
        d2l_ds2_k,
        max_ulps = 0,
        "scaling by a power of two is exact while every intermediate stays normal"
      );
    }
  }

  mod helpers {
    use ndarray::{Array1, array};
    use treetime_utils::cast::ToF64Exact;

    const JC69_DECAY_RATE: f64 = 4.0 / 3.0;

    pub(super) fn jc69_rate_spectrum() -> Array1<f64> {
      array![0.0, -JC69_DECAY_RATE, -JC69_DECAY_RATE, -JC69_DECAY_RATE]
    }

    pub(super) fn closed_form_in_sqrt_t(
      sqrt_t: f64,
      (stationary_weight, decaying_weight): (f64, f64),
      indel_count: usize,
      mu: f64,
    ) -> SqrtDerivatives {
      let beta = JC69_DECAY_RATE;
      let t = sqrt_t * sqrt_t;
      let decaying = decaying_weight * (-beta * t).exp();
      let site_lh = stationary_weight + decaying;
      let d_site_lh = -2.0 * beta * sqrt_t * decaying;
      let d2_site_lh = (4.0 * beta * beta * t - 2.0 * beta) * decaying;
      let substitution_first = d_site_lh / site_lh;
      let substitution_second = d2_site_lh / site_lh - substitution_first.powi(2);
      let indel_count = indel_count.to_f64_exact();
      let indel_first = 2.0 * indel_count / sqrt_t - 2.0 * mu * sqrt_t;
      let indel_second = -2.0 * indel_count / t - 2.0 * mu;
      SqrtDerivatives {
        first: substitution_first + indel_first,
        second: substitution_second + indel_second,
      }
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) struct SqrtDerivatives {
      pub(super) first: f64,
      pub(super) second: f64,
    }
  }

  mod generators {
    use proptest::prelude::*;

    pub(super) fn gen_s() -> impl Strategy<Value = f64> {
      1e-6_f64..1e3_f64
    }

    pub(super) fn gen_dl_dt() -> impl Strategy<Value = f64> {
      gen_signed(1e-3_f64..1e6_f64)
    }

    pub(super) fn gen_d2l_dt2() -> impl Strategy<Value = f64> {
      gen_signed(1e-3_f64..1e8_f64)
    }

    fn gen_signed(magnitude: impl Strategy<Value = f64>) -> impl Strategy<Value = f64> {
      (prop_oneof![Just(-1.0_f64), Just(1.0_f64)], magnitude).prop_map(|(sign, magnitude)| sign * magnitude)
    }

    pub(super) fn gen_same_sign_quadratic() -> impl Strategy<Value = (f64, f64)> {
      (
        prop_oneof![Just(-1.0_f64), Just(1.0_f64)],
        1e-3_f64..1e6_f64,
        1e-3_f64..1e8_f64,
      )
        .prop_map(|(sign, alpha, gamma)| (sign * alpha, sign * gamma))
    }

    pub(super) fn gen_power_of_two() -> impl Strategy<Value = f64> {
      (prop_oneof![Just(-1.0_f64), Just(1.0_f64)], -10_i32..=10_i32)
        .prop_map(|(sign, exponent)| sign * 2.0_f64.powi(exponent))
    }
  }
}
