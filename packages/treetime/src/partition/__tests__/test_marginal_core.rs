#[cfg(test)]
mod tests {
  use crate::partition::marginal_core::{
    forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_1d_inplace,
  };
  use ndarray::{Array1, array};
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::{pretty_assert_abs_diff_eq, prop_assert_ulps_eq};

  /// A matching degenerate scale belongs to the removed child factor and
  /// therefore cancels from the cavity message. See
  /// https://doi.org/10.1007/s00239-020-09982-w.
  #[test]
  fn test_marginal_core_forward_log_lh_remove_child_cancels_matching_neg_infinity() {
    let actual = forward_log_lh_remove_child(f64::NEG_INFINITY, f64::NEG_INFINITY);
    pretty_assert_abs_diff_eq!(0.0, actual, epsilon = 1e-12);
  }

  /// A degenerate normalization corresponds to the uniform fallback and has no
  /// finite scale contribution to a forward conditional message.
  #[test]
  fn test_marginal_core_forward_log_lh_add_normalization_ignores_neg_infinity() {
    let expected = -3.0;
    let actual = forward_log_lh_add_normalization(expected, f64::NEG_INFINITY);
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
  }

  /// Unexpected non-finite values remain visible instead of being treated as
  /// the negative-infinity fallback sentinel.
  #[test]
  fn test_marginal_core_forward_log_lh_add_normalization_propagates_nan() {
    let actual = forward_log_lh_add_normalization(-3.0, f64::NAN);
    assert!(actual.is_nan(), "NaN normalization must propagate, got {actual}");
  }

  /// Positive infinity is not the degenerate uniform-fallback sentinel.
  #[test]
  fn test_marginal_core_forward_log_lh_add_normalization_propagates_positive_infinity() {
    let actual = forward_log_lh_add_normalization(-3.0, f64::INFINITY);
    assert!(
      actual.is_infinite() && actual.is_sign_positive(),
      "expected positive infinity, got {actual}"
    );
  }

  /// A non-matching negative-infinity scale remains an impossible factor.
  #[test]
  fn test_marginal_core_forward_log_lh_remove_child_preserves_unmatched_neg_infinity() {
    let actual = forward_log_lh_remove_child(f64::NEG_INFINITY, -3.0);
    assert!(
      actual.is_infinite() && actual.is_sign_negative(),
      "expected negative infinity, got {actual}"
    );
  }

  /// Removing a negative-infinity scale from a finite aggregate remains
  /// observable as positive infinity rather than being mistaken for matching
  /// fallback sentinels.
  #[test]
  fn test_marginal_core_forward_log_lh_remove_child_preserves_unmatched_positive_infinity() {
    let actual = forward_log_lh_remove_child(-3.0, f64::NEG_INFINITY);
    assert!(
      actual.is_infinite() && actual.is_sign_positive(),
      "expected positive infinity, got {actual}"
    );
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    /// Finite forward normalization retains ordinary scale arithmetic.
    #[test]
    fn test_prop_marginal_core_forward_log_lh_matches_finite_arithmetic(
      node_log_lh in -1e6_f64..0.0,
      child_log_lh in -1e6_f64..0.0,
      normalization in -1e6_f64..0.0,
    ) {
      let cavity = forward_log_lh_remove_child(node_log_lh, child_log_lh);
      let actual = forward_log_lh_add_normalization(cavity, normalization);
      let expected = node_log_lh - child_log_lh + normalization;
      prop_assert_ulps_eq!(expected, actual, max_ulps = 0);
    }

    /// The degenerate normalization sentinel is the additive identity for
    /// forward conditional-message scales.
    #[test]
    fn test_prop_marginal_core_forward_log_lh_neg_infinity_is_neutral(log_lh in -1e6_f64..1e6_f64) {
      let actual = forward_log_lh_add_normalization(log_lh, f64::NEG_INFINITY);
      prop_assert_ulps_eq!(log_lh, actual, max_ulps = 0);
    }
  }

  // Valid (positive, finite) norm: distribution is normalized to sum 1 and the
  // contribution is weight * ln(norm).
  #[rustfmt::skip]
  #[rstest]
  #[case::unweighted(    array![1.0, 3.0],      1.0, array![0.25, 0.75],      4.0_f64.ln())]
  #[case::weighted(      array![1.0, 1.0],      5.0, array![0.5,  0.5],  5.0 * 2.0_f64.ln())]
  #[case::already_norm(  array![0.25, 0.75],    1.0, array![0.25, 0.75],      1.0_f64.ln())]
  #[case::three_states(  array![2.0, 2.0, 4.0], 1.0, array![0.25, 0.25, 0.5], 8.0_f64.ln())]
  #[trace]
  fn test_marginal_core_normalize_1d_inplace_valid(
    #[case] dis: Array1<f64>,
    #[case] weight: f64,
    #[case] expected_dis: Array1<f64>,
    #[case] expected_ll: f64,
  ) {
    let mut dis = dis;
    let ll = normalize_1d_inplace(&mut dis, weight);
    pretty_assert_abs_diff_eq!(expected_dis, dis, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(expected_ll, ll, epsilon = 1e-12);
  }

  // Non-positive or non-finite norm: fall back to a uniform distribution and
  // contribute NEG_INFINITY (unweighted, regardless of the weight argument).
  #[rustfmt::skip]
  #[rstest]
  #[case::zero_norm(     array![0.0, 0.0, 0.0],          3.0)]
  #[case::infinite_norm( array![f64::INFINITY, 1.0],     2.0)]
  #[case::nan_norm(      array![f64::NAN, 1.0],          1.0)]
  #[trace]
  fn test_marginal_core_normalize_1d_inplace_fallback(#[case] dis: Array1<f64>, #[case] weight: f64) {
    let mut dis = dis;
    let n = dis.len();
    let ll = normalize_1d_inplace(&mut dis, weight);
    let expected_uniform = Array1::from_elem(n, 1.0 / n as f64);
    pretty_assert_abs_diff_eq!(expected_uniform, dis, epsilon = 1e-12);
    assert!(ll.is_infinite() && ll.is_sign_negative(), "expected NEG_INFINITY, got {ll}");
  }
}
