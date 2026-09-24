#![allow(
  clippy::as_conversions,
  reason = "small index and state counts convert to f64 exactly"
)]

#[cfg(test)]
mod tests {
  use crate::array::softmax_with_log_norm::*;
  use crate::pretty_assert_neg_inf;
  use approx::assert_ulps_eq;
  use ndarray::Array1;
  use ndarray::array;
  use rstest::rstest;

  const NEG_INF: f64 = f64::NEG_INFINITY;

  #[rustfmt::skip]
  #[rstest]
  #[case::uniform_4(             Array1::from_elem(4, 0.25_f64.ln()))]
  #[case::descending(            array![0.0, -1.0, -2.0, -3.0])]
  #[case::ascending(             array![-3.0, -2.0, -1.0, 0.0])]
  #[case::all_same_nonzero(      array![3.0, 3.0, 3.0, 3.0])]
  #[case::close_values(          array![-0.001, -0.002, -0.003, -0.004])]
  #[case::wide_spread(           array![0.0, -10.0, -20.0, -30.0])]
  #[case::dominant_rest_tiny(    array![0.0, -100.0, -100.0, -100.0])]
  #[case::single_state(          array![5.0])]
  #[case::two_equal(             array![0.0, 0.0])]
  #[case::two_asymmetric(        array![0.0, -5.0])]
  #[case::three_ascending(       array![1.0, 2.0, 3.0])]
  #[case::twenty_uniform(        Array1::from_elem(20, 0.0))]
  #[case::nuc_realistic_profile( array![-0.1, -2.3, -4.5, -6.7])]
  #[trace]
  fn test_softmax_with_log_norm_finite(#[case] input: Array1<f64>) {
    let (normalized, log_norm) = softmax_with_log_norm(input.view());

    helpers::assert_valid_distribution(&normalized);
    assert!(log_norm.is_finite());
    assert!(
      normalized.iter().all(|&p| p > 0.0),
      "all-finite input must have all-positive output"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::mixed_finite_neg_inf(    array![0.0, NEG_INF, -1.0, NEG_INF])]
  #[case::single_finite_rest_inf(  array![NEG_INF, -3.0, NEG_INF, NEG_INF])]
  #[case::first_finite_rest_inf(   array![0.0, NEG_INF, NEG_INF, NEG_INF])]
  #[case::last_finite_rest_inf(    array![NEG_INF, NEG_INF, NEG_INF, 7.0])]
  #[case::two_finite_two_inf(      array![-1.0, NEG_INF, -2.0, NEG_INF])]
  #[trace]
  fn test_softmax_with_log_norm_with_neg_inf(#[case] input: Array1<f64>) {
    let (normalized, log_norm) = softmax_with_log_norm(input.view());

    helpers::assert_valid_distribution(&normalized);
    assert!(log_norm.is_finite());

    let expected_zeros = input.mapv(|v| v == NEG_INF);
    let actual_zeros = normalized.mapv(|p| p == 0.0);
    assert_eq!(expected_zeros, actual_zeros, "zero positions must match -inf positions");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::n1(  1)]
  #[case::n2(  2)]
  #[case::n3(  3)]
  #[case::n4(  4)]
  #[case::n20( 20)]
  #[trace]
  fn test_softmax_with_log_norm_degenerate(#[case] n_states: usize) {
    let log_vec = Array1::from_elem(n_states, NEG_INF);
    let (normalized, log_norm) = softmax_with_log_norm(log_vec.view());

    let expected_uniform = Array1::from_elem(n_states, 1.0 / n_states as f64);
    assert_ulps_eq!(normalized, expected_uniform, max_ulps = 0);
    pretty_assert_neg_inf!(log_norm, "expected NEG_INFINITY, got {log_norm}");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::n2(  2)]
  #[case::n4(  4)]
  #[case::n20( 20)]
  #[trace]
  fn test_softmax_with_log_norm_uniform_input(#[case] n_states: usize) {
    let log_vec = Array1::from_elem(n_states, 0.0);
    let (normalized, log_norm) = softmax_with_log_norm(log_vec.view());

    let expected_uniform = Array1::from_elem(n_states, 1.0 / n_states as f64);
    assert_ulps_eq!(normalized, expected_uniform, max_ulps = 0);
    assert!(log_norm.is_finite());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::shift_neg_1000( -1000.0)]
  #[case::shift_neg_100(   -100.0)]
  #[case::shift_pos_100(    100.0)]
  #[case::shift_pos_1000(  1000.0)]
  #[trace]
  fn test_softmax_with_log_norm_shift_invariant(#[case] offset: f64) {
    let base = array![0.0, -1.0, -2.0, -3.0];
    let shifted = base.mapv(|v| v + offset);

    let (base_norm, base_log_norm) = softmax_with_log_norm(base.view());
    let (shifted_norm, shifted_log_norm) = softmax_with_log_norm(shifted.view());

    helpers::assert_valid_distribution(&shifted_norm);
    assert_ulps_eq!(shifted_norm, base_norm, max_ulps = 2);
    assert_ulps_eq!(shifted_log_norm, base_log_norm + offset, max_ulps = 2);
  }

  mod stress {
    use super::*;

    #[test]
    fn test_softmax_with_log_norm_near_exp_overflow_boundary() {
      let input = array![709.0, 708.0, 707.0, 706.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm should be finite, got {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_near_exp_underflow_boundary() {
      let input = array![-744.0, -745.0, -746.0, -747.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm should be finite, got {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_at_f64_max_log() {
      let input = array![709.7, 709.6, 709.5, 709.4];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm overflowed: {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_log_norm_overflow() {
      let input = array![f64::MAX.ln() + 1.0, f64::MAX.ln(), f64::MAX.ln() - 1.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm={log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_extreme_spread() {
      let input = array![0.0, -1000.0, -1000.0, -1000.0];
      let (normalized, _log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert_eq!(normalized, array![1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_softmax_with_log_norm_alternating_extremes() {
      let input = array![500.0, -500.0, 500.0, -500.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite());
      assert_eq!(normalized, array![0.5, 0.0, 0.5, 0.0]);
    }

    #[test]
    fn test_softmax_with_log_norm_near_identical_at_high_magnitude() {
      let input = array![1e15, 1e15 + 1.0, 1e15 + 2.0, 1e15 + 3.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());

      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite());

      assert!(normalized[3] > normalized[2]);
      assert!(normalized[2] > normalized[1]);
      assert!(normalized[1] > normalized[0]);
    }

    #[test]
    fn test_softmax_with_log_norm_near_identical_catastrophic_cancellation() {
      let base = 1e15;
      let eps = base * f64::EPSILON;
      let input = array![base, base + eps, base + 2.0 * eps, base + 3.0 * eps];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite());
    }

    #[test]
    fn test_softmax_with_log_norm_nan_must_propagate() {
      let input = array![0.0, f64::NAN, -1.0, -2.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      assert!(
        normalized.iter().all(|v| v.is_nan()),
        "NaN input must propagate to every probability: {normalized}"
      );
      assert!(log_norm.is_nan(), "NaN input must propagate to log_norm: {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_all_nan_must_propagate() {
      let input = Array1::from_elem(4, f64::NAN);
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      assert!(
        normalized.iter().all(|v| v.is_nan()),
        "all-NaN input must propagate to every probability: {normalized}"
      );
      assert!(log_norm.is_nan(), "all-NaN input must propagate to log_norm: {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_positive_infinity_must_handle() {
      let input = array![0.0, f64::INFINITY, -1.0, -2.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert_eq!(normalized, array![0.0, 1.0, 0.0, 0.0]);
      assert!(log_norm.is_infinite() && log_norm.is_sign_positive());
    }
  }

  mod helpers {
    use approx::assert_abs_diff_eq;
    use ndarray::Array1;

    pub(super) fn assert_valid_distribution(normalized: &Array1<f64>) {
      assert!(
        normalized.iter().all(|&v| v >= 0.0),
        "all probabilities must be non-negative: {normalized}"
      );
      assert!(
        normalized.iter().all(|&v| v.is_finite()),
        "all probabilities must be finite: {normalized}"
      );
      assert_abs_diff_eq!(normalized.sum(), 1.0, epsilon = 1e-15);
    }
  }
}
