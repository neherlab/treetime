#[cfg(test)]
mod tests {
  use crate::optimize::method_newton::sqrt_step_lower_bound;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_optimize_method_step_clamping_sqrt_at_zero() {
    pretty_assert_abs_diff_eq!(sqrt_step_lower_bound(0.0), -1.0, epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::s_0_3( 0.3, 0.3 - (0.09_f64 + 1.0).sqrt())]
  #[case::s_1_0( 1.0, 1.0 - 2.0_f64.sqrt())]
  #[case::s_3_0( 3.0, 3.0 - 10.0_f64.sqrt())]
  #[case::s_10(  10.0, 10.0 - 101.0_f64.sqrt())]
  #[trace]
  fn test_optimize_method_step_clamping_sqrt_analytical(#[case] s: f64, #[case] expected: f64) {
    pretty_assert_abs_diff_eq!(sqrt_step_lower_bound(s), expected, epsilon = 1e-14);
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
    pretty_assert_abs_diff_eq!(t_new - t_old, 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_optimize_method_step_clamping_sqrt_large_s_asymptote() {
    let s = 1000.0;
    let bound = sqrt_step_lower_bound(s);
    let asymptote = -1.0 / (2.0 * s);
    pretty_assert_abs_diff_eq!(bound, asymptote, epsilon = 1e-9);
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

  mod generators {
    use proptest::prelude::*;
    pub(super) fn gen_s() -> impl Strategy<Value = f64> {
      0.0_f64..1e3_f64
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
  }
}
