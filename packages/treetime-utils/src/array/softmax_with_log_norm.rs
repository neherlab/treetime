#[cfg(test)]
mod __tests__;

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
