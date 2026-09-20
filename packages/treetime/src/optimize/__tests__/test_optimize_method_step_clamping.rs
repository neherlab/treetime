#[cfg(test)]
mod tests {
  use crate::optimize::method_newton::{log_step_lower_bound, sqrt_step_lower_bound};
  use approx::assert_abs_diff_eq;
  use rstest::rstest;

  #[test]
  fn test_optimize_method_step_clamping_sqrt_at_zero() {
    assert_abs_diff_eq!(sqrt_step_lower_bound(0.0), -1.0, epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::s_0_3( 0.3, 0.3 - (0.09_f64 + 1.0).sqrt())]
  #[case::s_1_0( 1.0, 1.0 - 2.0_f64.sqrt())]
  #[case::s_3_0( 3.0, 3.0 - 10.0_f64.sqrt())]
  #[case::s_10(  10.0, 10.0 - 101.0_f64.sqrt())]
  #[trace]
  fn test_optimize_method_step_clamping_sqrt_analytical(#[case] s: f64, #[case] expected: f64) {
    assert_abs_diff_eq!(sqrt_step_lower_bound(s), expected, epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(  0.0)]
  #[case::tiny(  1e-10)]
  #[case::small( 0.01)]
  #[case::medium(0.3)]
  #[case::one(   1.0)]
  #[case::large( 10.0)]
  #[case::huge(  1000.0)]
  #[trace]
  fn test_optimize_method_step_clamping_sqrt_always_negative(#[case] s: f64) {
    assert!(sqrt_step_lower_bound(s) < 0.0);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tiny(  1e-6)]
  #[case::small( 0.01)]
  #[case::medium(0.3)]
  #[case::one(   1.0)]
  #[case::large( 5.0)]
  #[case::huge(  100.0)]
  #[trace]
  fn test_optimize_method_step_clamping_sqrt_produces_unit_t_increase(#[case] s: f64) {
    let delta_s = sqrt_step_lower_bound(s);
    let s_new = s - delta_s;
    let t_old = s * s;
    let t_new = s_new * s_new;
    assert_abs_diff_eq!(t_new - t_old, 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_optimize_method_step_clamping_sqrt_large_s_asymptote() {
    let s = 1000.0;
    let bound = sqrt_step_lower_bound(s);
    let asymptote = -1.0 / (2.0 * s);
    assert_abs_diff_eq!(bound, asymptote, epsilon = 1e-9);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small( 0.01)]
  #[case::medium(0.3)]
  #[case::one(   1.0)]
  #[case::large( 10.0)]
  #[trace]
  fn test_optimize_method_step_clamping_sqrt_strictly_greater_than_minus_one(#[case] s: f64) {
    let bound = sqrt_step_lower_bound(s);
    assert!(
      bound > -1.0,
      "At s={s}, bound ({bound}) should be > -1.0"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_1( 0.1,  -(10.0_f64).ln_1p())]
  #[case::t_1_0( 1.0,  -(1.0_f64).ln_1p())]
  #[case::t_3_0( 3.0,  -((1.0_f64 / 3.0).ln_1p()))]
  #[case::t_10(  10.0, -(0.1_f64).ln_1p())]
  #[trace]
  fn test_optimize_method_step_clamping_log_analytical(#[case] t: f64, #[case] expected: f64) {
    assert_abs_diff_eq!(log_step_lower_bound(t), expected, epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tiny(  1e-10)]
  #[case::small( 0.01)]
  #[case::medium(0.3)]
  #[case::one(   1.0)]
  #[case::large( 10.0)]
  #[case::huge(  1000.0)]
  #[trace]
  fn test_optimize_method_step_clamping_log_always_negative(#[case] t: f64) {
    assert!(log_step_lower_bound(t) < 0.0);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tiny(  1e-6)]
  #[case::small( 0.01)]
  #[case::medium(0.3)]
  #[case::one(   1.0)]
  #[case::large( 5.0)]
  #[case::huge(  100.0)]
  #[trace]
  fn test_optimize_method_step_clamping_log_produces_unit_t_increase(#[case] t: f64) {
    let delta_u = log_step_lower_bound(t);
    let u_new = t.ln() - delta_u;
    let t_new = u_new.exp();
    assert_abs_diff_eq!(t_new - t, 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_optimize_method_step_clamping_log_large_t_asymptote() {
    let t = 1000.0;
    let bound = log_step_lower_bound(t);
    let asymptote = -1.0 / t;
    assert_abs_diff_eq!(bound, asymptote, epsilon = 1e-6);
  }

  mod generators {
    use proptest::prelude::*;
    pub fn gen_s() -> impl Strategy<Value = f64> {
      0.0_f64..1e3_f64
    }
    pub fn gen_t() -> impl Strategy<Value = f64> {
      1e-10_f64..1e6_f64
    }
  }

  use proptest::prelude::*;

  proptest! {
    #[test]
    fn test_prop_optimize_method_step_clamping_sqrt_non_positive(s in generators::gen_s()) {
      let bound = sqrt_step_lower_bound(s);
      prop_assert!(bound <= 0.0, "sqrt_step_lower_bound({s}) = {bound} should be <= 0");
    }

    #[test]
    fn test_prop_optimize_method_step_clamping_sqrt_unit_t_increase(s in generators::gen_s()) {
      let delta_s = sqrt_step_lower_bound(s);
      let s_new = s - delta_s;
      let delta_t = s_new * s_new - s * s;
      prop_assert!(
        (delta_t - 1.0).abs() < 1e-6,
        "sqrt step from s={s} should yield delta_t = 1.0, got {delta_t}"
      );
    }

    #[test]
    fn test_prop_optimize_method_step_clamping_log_negative(t in generators::gen_t()) {
      let bound = log_step_lower_bound(t);
      prop_assert!(bound < 0.0, "log_step_lower_bound({t}) = {bound} should be < 0");
    }

    #[test]
    fn test_prop_optimize_method_step_clamping_log_unit_t_increase(t in generators::gen_t()) {
      let delta_u = log_step_lower_bound(t);
      let t_new = (t.ln() - delta_u).exp();
      let delta_t = t_new - t;
      prop_assert!(
        (delta_t - 1.0).abs() < 1e-6,
        "log step from t={t} should yield delta_t = 1.0, got {delta_t}"
      );
    }
  }
}
