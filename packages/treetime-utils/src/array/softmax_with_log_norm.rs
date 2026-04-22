use ndarray::{Array1, ArrayView1};
use ndarray_stats::QuantileExt;

/// Fused softmax + log-sum-exp (LSE): returns `(softmax(x), LSE(x))` from log-space inputs.
///
/// # Softmax
///
/// The softmax function converts a vector of real numbers into a probability distribution:
///
/// $$\text{softmax}(\mathbf{x})_i = \frac{e^{x_i}}{\sum_j e^{x_j}}$$
///
/// When inputs are log-probabilities (or log-likelihoods), softmax exponentiates and
/// normalizes them into proper probabilities that sum to 1.
///
/// # Log-Sum-Exp (LSE)
///
/// The log-sum-exp function computes:
///
/// $$\text{LSE}(\mathbf{x}) = \ln\left(\sum_i e^{x_i}\right)$$
///
/// For log-probabilities, LSE equals the log of the normalization constant (log partition
/// function). This is the value subtracted to normalize:
/// $\text{softmax}(\mathbf{x})_i = e^{x_i - \text{LSE}(\mathbf{x})}$.
///
/// # Numerical Stability
///
/// Naive computation of $e^x$ overflows for $x > 709.78$ (`f64::MAX.ln()`) and underflows
/// to zero for $x < -745.13$ (`f64::MIN_POSITIVE.ln()`). Log-probabilities often span this
/// range (e.g., summing many small probabilities produces deeply negative log values).
///
/// The standard trick subtracts the maximum before exponentiation:
///
/// $$\text{LSE}(\mathbf{x}) = \max(\mathbf{x}) + \ln\left(\sum_i e^{x_i - \max(\mathbf{x})}\right)$$
///
/// After shifting by $\max(\mathbf{x})$, the largest exponent is $e^0 = 1$, keeping all
/// values in safe range. The shift does not change the result: $\ln(e^c \cdot X) = c + \ln(X)$.
///
/// # Fusion
///
/// Computing softmax and LSE separately would duplicate work. The decomposed softmax
/// $e^{x_i - \text{LSE}(\mathbf{x})}$ requires computing $e^{x_i - \max}$ twice and
/// introduces an extra $\ln \to \exp$ roundtrip.
///
/// Fused computation shares the intermediate $e^{x_i - \max}$:
///
/// $$s_i = e^{x_i - \max}, \quad S = \sum_i s_i$$
/// $$\text{LSE}(\mathbf{x}) = \max + \ln(S), \quad \text{softmax}(\mathbf{x})_i = \frac{s_i}{S}$$
///
/// Direct division $s_i / S$ preserves tighter sum-to-one invariants than subtracting
/// a recomputed log normalization.
///
/// # Returns
///
/// `(softmax, log_norm)` where:
/// - `softmax`: normalized probability distribution (sums to 1.0)
/// - `log_norm`: log-sum-exp value (log of normalization constant)
///
/// # Edge Cases
///
/// - **Empty input**: returns `([], -inf)`. Zero-element sum has log zero = -inf.
/// - **NaN present**: returns `(NaN, NaN)`. Per IEEE 754, NaN propagates.
/// - **+inf present**: one-hot on infinite positions, `log_norm = +inf`. Infinite
///   log-probability dominates; mass distributed uniformly across all +inf positions.
/// - **All -inf**: uniform distribution, `log_norm = -inf`. Zero probability everywhere,
///   but softmax must be a valid distribution, so uniform is the only symmetric choice.
///
/// # Example
///
/// ```
/// use ndarray::array;
/// use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm;
///
/// let log_probs = array![-1.0, -2.0, -3.0];
/// let (probs, log_norm) = softmax_with_log_norm(log_probs.view());
///
/// // probs sums to 1.0, highest probability at index 0 (highest log-prob)
/// assert!((probs.sum() - 1.0).abs() < 1e-15);
/// assert!(probs[0] > probs[1] && probs[1] > probs[2]);
/// ```
pub fn softmax_with_log_norm(log_vec: ArrayView1<'_, f64>) -> (Array1<f64>, f64) {
  let n = log_vec.len();

  // LSE shift constant: subtract max before exp to prevent overflow.
  // QuantileExt::max() errors on empty or NaN-containing arrays.
  let max_val = match log_vec.max() {
    Ok(&max) => max,
    Err(_) if n == 0 => return (Array1::zeros(0), f64::NEG_INFINITY),
    Err(_) => return (Array1::from_elem(n, f64::NAN), f64::NAN),
  };

  // +inf dominates: uniform weight across all infinite positions
  if max_val == f64::INFINITY {
    let inf_count = log_vec.iter().filter(|&&v| v == f64::INFINITY).count();
    let prob = 1.0 / inf_count as f64;
    return (
      log_vec.mapv(|v| if v == f64::INFINITY { prob } else { 0.0 }),
      f64::INFINITY,
    );
  }

  // All -inf: uniform fallback
  if !max_val.is_finite() {
    return (Array1::from_elem(n, 1.0 / n as f64), f64::NEG_INFINITY);
  }

  // Shared intermediate: exp(x - max) used by both LSE and softmax
  let shifted = log_vec.mapv(|v| (v - max_val).exp());
  let shifted_sum = shifted.sum();

  // LSE(x) = max + ln(sum(exp(x - max)))
  let log_norm = max_val + shifted_sum.ln();

  // softmax(x) = exp(x - max) / sum(exp(x - max))
  let normalized = shifted / shifted_sum;

  (normalized, log_norm)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_neg_inf;
  use approx::assert_ulps_eq;
  use ndarray::array;
  use rstest::rstest;

  const NEG_INF: f64 = f64::NEG_INFINITY;

  /// All-finite inputs: softmax output has no zeros and all probabilities are positive.
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

  /// Inputs containing -inf: positions with -inf get probability 0.0, others share the mass.
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

    // Zeros must appear exactly where input is -inf
    let expected_zeros = input.mapv(|v| v == NEG_INF);
    let actual_zeros = normalized.mapv(|p| p == 0.0);
    assert_eq!(expected_zeros, actual_zeros, "zero positions must match -inf positions");
  }

  /// All-inf degenerate input: uniform fallback, log_norm = -inf.
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

  /// Equal finite inputs produce exact uniform output.
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

  /// Shift invariance: adding a constant to all inputs preserves softmax,
  /// shifts log_norm by the same constant.
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

  /// Numerical stress tests: attempt to break LSE stability.
  mod stress {
    use super::*;

    #[test]
    fn test_softmax_with_log_norm_near_exp_overflow_boundary() {
      // exp(709.78) ~ f64::MAX, exp(710) overflows.
      // After shift by max=709.0, exponents are [0, -1, -2, -3]. Safe.
      let input = array![709.0, 708.0, 707.0, 706.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm should be finite, got {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_near_exp_underflow_boundary() {
      // exp(-745) underflows to 0. After shift by max=-744.0,
      // exponents are [0, -1, -2, -3]. Safe.
      let input = array![-744.0, -745.0, -746.0, -747.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm should be finite, got {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_at_f64_max_log() {
      // log_norm = max + ln(sum), if max is near ln(f64::MAX) ~ 709.78,
      // then log_norm could overflow. Push to the edge.
      let input = array![709.7, 709.6, 709.5, 709.4];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite(), "log_norm overflowed: {log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_log_norm_overflow() {
      // log_norm = max + ln(sum(exp(x - max))). When max > f64::MAX.ln() ~ 709.78
      // and sum > 1, log_norm can exceed f64::MAX.ln(). But log_norm is a log-value,
      // not an exponentiated value, so it remains finite well beyond exp's overflow.
      // Push to actual f64 overflow: max near f64::MAX itself.
      let input = array![f64::MAX.ln() + 1.0, f64::MAX.ln(), f64::MAX.ln() - 1.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      // log_norm ~ 710.78 + ln(1 + e^-1 + e^-2) ~ 711.23, still finite.
      assert!(log_norm.is_finite(), "log_norm={log_norm}");
    }

    #[test]
    fn test_softmax_with_log_norm_extreme_spread() {
      // One dominant state at 0, rest at -1000. The non-dominant states
      // underflow to exp(-1000) = 0 after shift, producing exact one-hot.
      let input = array![0.0, -1000.0, -1000.0, -1000.0];
      let (normalized, _log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert_eq!(normalized, array![1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_softmax_with_log_norm_alternating_extremes() {
      // Alternating high/low. After shift by 500, exponents are [0, -1000, 0, -1000].
      // Low positions underflow to 0, high positions share mass equally.
      let input = array![500.0, -500.0, 500.0, -500.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite());
      assert_eq!(normalized, array![0.5, 0.0, 0.5, 0.0]);
    }

    #[test]
    fn test_softmax_with_log_norm_near_identical_at_high_magnitude() {
      // Values differ by 1.0 at 1e15 magnitude. Machine epsilon at 1e15 is ~0.22,
      // so differences of 1.0 are ~4.5 ULPs apart. Precision should be preserved.
      let input = array![1e15, 1e15 + 1.0, 1e15 + 2.0, 1e15 + 3.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());

      helpers::assert_valid_distribution(&normalized);
      assert!(log_norm.is_finite());

      // Monotonicity: higher input -> higher probability
      assert!(normalized[3] > normalized[2]);
      assert!(normalized[2] > normalized[1]);
      assert!(normalized[1] > normalized[0]);
    }

    #[test]
    fn test_softmax_with_log_norm_near_identical_catastrophic_cancellation() {
      // Values differ by less than machine epsilon * magnitude.
      // x - max loses all significant digits, producing junk relative
      // probabilities. The output is still a valid distribution (sums
      // to 1, non-negative) but does not reflect true relative weights.
      // This is inherent to finite-precision arithmetic.
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
      let has_nan = normalized.iter().any(|v| v.is_nan()) || log_norm.is_nan();
      assert!(
        has_nan,
        "NaN input must propagate: normalized={normalized}, log_norm={log_norm}"
      );
    }

    #[test]
    fn test_softmax_with_log_norm_all_nan_must_propagate() {
      let input = Array1::from_elem(4, f64::NAN);
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      let has_nan = normalized.iter().any(|v| v.is_nan()) || log_norm.is_nan();
      assert!(
        has_nan,
        "all-NaN input must propagate: normalized={normalized}, log_norm={log_norm}"
      );
    }

    #[test]
    fn test_softmax_with_log_norm_positive_infinity_must_handle() {
      // +inf log-probability means infinite dominance for that state.
      // Correct output: exact one-hot on the +inf state, log_norm = +inf.
      let input = array![0.0, f64::INFINITY, -1.0, -2.0];
      let (normalized, log_norm) = softmax_with_log_norm(input.view());
      helpers::assert_valid_distribution(&normalized);
      assert_eq!(normalized, array![0.0, 1.0, 0.0, 0.0]);
      assert!(log_norm == f64::INFINITY);
    }
  }

  mod helpers {
    use approx::assert_abs_diff_eq;
    use ndarray::Array1;

    pub fn assert_valid_distribution(normalized: &Array1<f64>) {
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
