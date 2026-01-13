use crate::traits::MultiplyAlgo;
use ndarray::Array1;

/// Threshold below which we normalize to prevent underflow.
/// This is much smaller than 1.0 to avoid unnecessary normalizations.
const UNDERFLOW_THRESHOLD: f64 = 1e-100;

fn array_max(arr: &Array1<f64>) -> f64 {
  arr.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn normalize_and_extract_scale(arr: &Array1<f64>) -> (Array1<f64>, f64) {
  let max_val = array_max(arr);
  let log_scale = if max_val > 0.0 && max_val.is_finite() {
    max_val.ln()
  } else {
    f64::NEG_INFINITY
  };
  let normalized = if max_val > 0.0 && max_val.is_finite() {
    arr.mapv(|v| v / max_val)
  } else {
    arr.clone()
  };
  (normalized, log_scale)
}

/// Multiply multiple distributions with lazy normalization.
///
/// Returns (normalized_shape, log_scale) where the full result is:
/// `normalized_shape * exp(log_scale)`
///
/// Invariant: max(normalized_shape) = 1.0
///
/// Unlike aggressive normalization (which normalizes after every step),
/// this function only normalizes when the maximum value drops below
/// `UNDERFLOW_THRESHOLD` to prevent numerical underflow. This reduces
/// accumulated rounding errors and improves both performance and accuracy.
pub fn multiply_many_lazy_normalize(distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
  if distributions.is_empty() {
    return (Array1::zeros(0), f64::NEG_INFINITY);
  }

  if distributions.len() == 1 {
    return normalize_and_extract_scale(distributions[0]);
  }

  let n = distributions[0].len();
  let mut product = distributions[0].clone();
  let mut accumulated_log_scale = 0.0;

  // Initial check
  let max_val = array_max(&product);
  if max_val <= 0.0 || !max_val.is_finite() {
    return (Array1::zeros(n), f64::NEG_INFINITY);
  }

  for dist in distributions.iter().skip(1) {
    product = &product * *dist;

    let max_val = array_max(&product);
    if max_val <= 0.0 || !max_val.is_finite() {
      return (Array1::zeros(n), f64::NEG_INFINITY);
    }

    // Only normalize if approaching underflow
    if max_val < UNDERFLOW_THRESHOLD {
      accumulated_log_scale += max_val.ln();
      product.mapv_inplace(|v| v / max_val);
    }
  }

  // Final normalization to ensure invariant max = 1.0
  let max_val = array_max(&product);
  if max_val <= 0.0 || !max_val.is_finite() {
    return (Array1::zeros(n), f64::NEG_INFINITY);
  }
  accumulated_log_scale += max_val.ln();
  product.mapv_inplace(|v| v / max_val);

  (product, accumulated_log_scale)
}

/// Multiply multiple distributions with aggressive normalization (legacy behavior).
///
/// Returns (normalized_shape, log_scale) where the full result is:
/// `normalized_shape * exp(log_scale)`
///
/// Invariant: max(normalized_shape) = 1.0
///
/// This function normalizes after EVERY pairwise multiplication.
/// Prefer `multiply_many_lazy_normalize` for better accuracy.
pub fn multiply_many(distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
  if distributions.is_empty() {
    return (Array1::zeros(0), f64::NEG_INFINITY);
  }

  if distributions.len() == 1 {
    return normalize_and_extract_scale(distributions[0]);
  }

  let n = distributions[0].len();
  let mut accumulated_log_scale = 0.0;
  let mut normalized_result = distributions[0].clone();

  let max_val = array_max(&normalized_result);
  if max_val <= 0.0 || !max_val.is_finite() {
    return (Array1::zeros(n), f64::NEG_INFINITY);
  }
  accumulated_log_scale += max_val.ln();
  normalized_result.mapv_inplace(|v| v / max_val);

  for dist in distributions.iter().skip(1) {
    normalized_result = &normalized_result * *dist;

    let max_val = array_max(&normalized_result);
    if max_val <= 0.0 || !max_val.is_finite() {
      return (Array1::zeros(n), f64::NEG_INFINITY);
    }

    accumulated_log_scale += max_val.ln();
    normalized_result.mapv_inplace(|v| v / max_val);
  }

  (normalized_result, accumulated_log_scale)
}

/// Naive multiplication without log-scale tracking.
///
/// Warning: Can underflow to zero for many small distributions.
/// Use `multiply_many_lazy_normalize` for underflow-resistant multiplication.
pub fn multiply_many_naive(distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
  if distributions.is_empty() {
    return (Array1::zeros(0), f64::NEG_INFINITY);
  }

  let n = distributions[0].len();
  let mut result = Array1::ones(n);

  for dist in distributions {
    result = &result * *dist;
  }

  let max_val = array_max(&result);
  let log_scale = if max_val > 0.0 && max_val.is_finite() {
    max_val.ln()
  } else {
    f64::NEG_INFINITY
  };

  if max_val > 0.0 && max_val.is_finite() {
    result.mapv_inplace(|v| v / max_val);
  }

  (result, log_scale)
}

/// Naive point-wise multiplication algorithm.
pub struct PointwiseMultiply;

impl MultiplyAlgo for PointwiseMultiply {
  fn name(&self) -> &'static str {
    "naive-multiplication"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
    multiply_many_naive(distributions)
  }
}

/// Log-scale aware multiplication with lazy normalization.
pub struct LogScaleMultiply;

impl MultiplyAlgo for LogScaleMultiply {
  fn name(&self) -> &'static str {
    "log-scale-multiplication"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
    multiply_many_lazy_normalize(distributions)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_relative_eq;
  use ndarray::array;

  #[test]
  fn test_multiply_many_empty() {
    let (shape, log_scale) = multiply_many_lazy_normalize(&[]);
    assert_eq!(shape.len(), 0);
    assert!(log_scale.is_infinite() && log_scale.is_sign_negative());
  }

  #[test]
  fn test_multiply_many_single() {
    let a = array![1.0, 4.0, 2.0];
    let (shape, log_scale) = multiply_many_lazy_normalize(&[&a]);

    assert_relative_eq!(shape[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(shape[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(shape[2], 0.5, epsilon = 1e-10);
    assert_relative_eq!(log_scale, 4.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_multiply_many_two_identical() {
    let a = array![1.0, 2.0, 1.0];
    let (shape, log_scale) = multiply_many_lazy_normalize(&[&a, &a]);

    assert_relative_eq!(shape[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(shape[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(shape[2], 0.25, epsilon = 1e-10);
    assert!(log_scale.is_finite());
  }

  #[test]
  fn test_multiply_many_underflow_resistance() {
    let n = 100;
    let small_val = 0.0001;
    let a = array![small_val];

    let refs: Vec<&Array1<f64>> = std::iter::repeat_n(&a, n).collect();
    let (shape, log_scale) = multiply_many_lazy_normalize(&refs);

    assert!(!shape.is_empty());
    assert!(log_scale.is_finite());

    let expected_log_scale = (n as f64) * small_val.ln();
    assert_relative_eq!(log_scale, expected_log_scale, epsilon = 1e-10);
  }

  #[test]
  fn test_multiply_many_handles_zero_max() {
    let a = array![0.0, 0.0, 0.0];
    let (shape, log_scale) = multiply_many_lazy_normalize(&[&a]);

    assert!(log_scale.is_infinite() && log_scale.is_sign_negative());
    assert_eq!(shape, a);
  }

  #[test]
  fn test_lazy_vs_aggressive_same_result_for_few_factors() {
    // For a small number of factors, both should give very similar results
    let a = array![1.0, 2.0, 1.0];
    let b = array![0.5, 1.0, 0.5];

    let (shape_lazy, log_lazy) = multiply_many_lazy_normalize(&[&a, &b]);
    let (shape_agg, log_agg) = multiply_many(&[&a, &b]);

    // Results should be very close
    for i in 0..shape_lazy.len() {
      assert_relative_eq!(shape_lazy[i], shape_agg[i], epsilon = 1e-10);
    }
    assert_relative_eq!(log_lazy, log_agg, epsilon = 1e-10);
  }

  #[test]
  fn test_lazy_normalize_fewer_operations() {
    // With large values that don't risk underflow, lazy should not normalize mid-computation
    // This is hard to test directly, but we can verify correctness
    let a = array![0.5, 1.0, 0.5];
    let b = array![0.5, 1.0, 0.5];
    let c = array![0.5, 1.0, 0.5];

    let (shape, log_scale) = multiply_many_lazy_normalize(&[&a, &b, &c]);

    // a * b * c at center = 1.0 * 1.0 * 1.0 = 1.0
    // a * b * c at edges = 0.5 * 0.5 * 0.5 = 0.125
    // max = 1.0, so normalized shape[1] = 1.0, shape[0] = shape[2] = 0.125
    assert_relative_eq!(shape[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(shape[0], 0.125, epsilon = 1e-10);
    assert_relative_eq!(shape[2], 0.125, epsilon = 1e-10);
    assert_relative_eq!(log_scale, 0.0, epsilon = 1e-10); // max was 1.0
  }
}
