#[cfg(test)]
mod tests {
  use crate::testing::assert::is_neg_inf;
  use crate::{pretty_assert_map_ulps_eq, pretty_assert_neg_inf, pretty_assert_relative_eq, pretty_assert_ulps_eq};
  use helpers::catch_panic_silent;
  use maplit::btreemap;
  use ndarray::array;
  use rstest::rstest;
  use std::panic::catch_unwind;

  #[test]
  fn test_is_neg_inf_accepts_only_negative_infinity() {
    assert!(is_neg_inf(f64::NEG_INFINITY));
    assert!(!is_neg_inf(f64::INFINITY));
    assert!(!is_neg_inf(f64::NAN));
    assert!(!is_neg_inf(-1.0));
    assert!(!is_neg_inf(-0.0));
    assert!(!is_neg_inf(0.0));
  }

  #[test]
  fn test_pretty_assert_neg_inf_accepts_negative_infinity() {
    let result = catch_unwind(|| pretty_assert_neg_inf!(f64::NEG_INFINITY));
    result.unwrap();
  }

  #[test]
  fn test_pretty_assert_neg_inf_rejects_positive_infinity() {
    let result = catch_panic_silent(|| pretty_assert_neg_inf!(f64::INFINITY));
    assert!(result.is_err());
  }

  #[test]
  fn test_pretty_assert_neg_inf_rejects_nan() {
    let result = catch_panic_silent(|| pretty_assert_neg_inf!(f64::NAN));
    assert!(result.is_err());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tiny_values_twice_apart(  1e-20,   2e-20)]
  #[case::zero_and_tiny_value(      0.0,     1e-300)]
  #[case::tiny_values_of_both_signs(-1e-300, 1e-300)]
  #[trace]
  fn test_pretty_assert_ulps_eq_rejects_tiny_values_more_than_max_ulps_apart(
    #[case] lhs: f64,
    #[case] rhs: f64,
  ) {
    let result = catch_panic_silent(|| pretty_assert_ulps_eq!(lhs, rhs, max_ulps = 4));
    assert!(result.is_err());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_and_negative_zero(0.0,   -0.0)]
  #[case::adjacent_tiny_values(  1e-20, 1e-20_f64.next_up())]
  #[case::adjacent_values(       1.0,   1.0_f64.next_up())]
  #[trace]
  fn test_pretty_assert_ulps_eq_accepts_values_within_max_ulps(#[case] lhs: f64, #[case] rhs: f64) {
    pretty_assert_ulps_eq!(lhs, rhs, max_ulps = 1);
  }

  #[test]
  fn test_pretty_assert_ulps_eq_accepts_values_within_the_given_epsilon() {
    pretty_assert_ulps_eq!(0.0, 1e-17, epsilon = 1e-16, max_ulps = 4);
  }

  #[test]
  fn test_pretty_assert_ulps_eq_rejects_arrays_of_tiny_values_far_apart() {
    let result = catch_panic_silent(|| pretty_assert_ulps_eq!(array![1.0, 1e-20], array![1.0, 2e-20], max_ulps = 4));
    assert!(result.is_err());
  }

  #[test]
  fn test_pretty_assert_map_ulps_eq_rejects_tiny_values_far_apart() {
    let result = catch_panic_silent(|| pretty_assert_map_ulps_eq!(btreemap! {1 => 1e-20}, btreemap! {1 => 2e-20}));
    assert!(result.is_err());
  }

  #[test]
  fn test_pretty_assert_relative_eq_rejects_tiny_values_far_apart() {
    let result = catch_panic_silent(|| pretty_assert_relative_eq!(1e-20, 2e-20, max_relative = 1e-10));
    assert!(result.is_err());
  }

  #[test]
  fn test_pretty_assert_relative_eq_accepts_tiny_values_within_max_relative() {
    pretty_assert_relative_eq!(1e-20, 1.000_000_000_01e-20, max_relative = 1e-10);
  }

  mod helpers {
    use parking_lot::Mutex;
    use std::panic::{AssertUnwindSafe, catch_unwind, set_hook, take_hook};
    use std::sync::LazyLock;
    use std::thread;

    static PANIC_HOOK_LOCK: LazyLock<Mutex<()>> = LazyLock::new(|| Mutex::new(()));

    pub(super) fn catch_panic_silent(f: impl FnOnce()) -> thread::Result<()> {
      let _guard = PANIC_HOOK_LOCK.lock();
      let original_hook = take_hook();
      set_hook(Box::new(|_| {}));
      let result = catch_unwind(AssertUnwindSafe(f));
      set_hook(original_hook);
      result
    }
  }
}
