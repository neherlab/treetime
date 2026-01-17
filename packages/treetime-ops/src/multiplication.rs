use crate::ScaledArray;
use crate::traits::MultiplyAlgo;
use ndarray::Array1;

/// Threshold below which we normalize to prevent underflow.
/// This is much smaller than 1.0 to avoid unnecessary normalizations.
const UNDERFLOW_THRESHOLD: f64 = 1e-100;

fn array_max(arr: &Array1<f64>) -> f64 {
  arr.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn normalize_and_extract_scale(arr: &Array1<f64>) -> ScaledArray {
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
  ScaledArray::new(normalized, log_scale)
}

/// Multiply multiple distributions with lazy normalization.
///
/// Returns `ScaledArray { normalized, log_scale }` where the full result is:
/// `normalized * exp(log_scale)`
///
/// Invariant: max(normalized) = 1.0
///
/// Unlike aggressive normalization (which normalizes after every step),
/// this function only normalizes when the maximum value drops below
/// `UNDERFLOW_THRESHOLD` to prevent numerical underflow. This reduces
/// accumulated rounding errors and improves both performance and accuracy.
pub fn multiply_many_lazy_normalize(distributions: &[&Array1<f64>]) -> ScaledArray {
  if distributions.is_empty() {
    return ScaledArray::empty(0);
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
    return ScaledArray::empty(n);
  }

  for dist in distributions.iter().skip(1) {
    product = &product * *dist;

    let max_val = array_max(&product);
    if max_val <= 0.0 || !max_val.is_finite() {
      return ScaledArray::empty(n);
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
    return ScaledArray::empty(n);
  }
  accumulated_log_scale += max_val.ln();
  product.mapv_inplace(|v| v / max_val);

  ScaledArray::new(product, accumulated_log_scale)
}

/// Multiply multiple distributions with aggressive normalization (legacy behavior).
///
/// Returns `ScaledArray { normalized, log_scale }` where the full result is:
/// `normalized * exp(log_scale)`
///
/// Invariant: max(normalized) = 1.0
///
/// This function normalizes after EVERY pairwise multiplication.
/// Prefer `multiply_many_lazy_normalize` for better accuracy.
pub fn multiply_many(distributions: &[&Array1<f64>]) -> ScaledArray {
  if distributions.is_empty() {
    return ScaledArray::empty(0);
  }

  if distributions.len() == 1 {
    return normalize_and_extract_scale(distributions[0]);
  }

  let n = distributions[0].len();
  let mut accumulated_log_scale = 0.0;
  let mut normalized_result = distributions[0].clone();

  let max_val = array_max(&normalized_result);
  if max_val <= 0.0 || !max_val.is_finite() {
    return ScaledArray::empty(n);
  }
  accumulated_log_scale += max_val.ln();
  normalized_result.mapv_inplace(|v| v / max_val);

  for dist in distributions.iter().skip(1) {
    normalized_result = &normalized_result * *dist;

    let max_val = array_max(&normalized_result);
    if max_val <= 0.0 || !max_val.is_finite() {
      return ScaledArray::empty(n);
    }

    accumulated_log_scale += max_val.ln();
    normalized_result.mapv_inplace(|v| v / max_val);
  }

  ScaledArray::new(normalized_result, accumulated_log_scale)
}

/// Naive multiplication without log-scale tracking.
///
/// Warning: Can underflow to zero for many small distributions.
/// Use `multiply_many_lazy_normalize` for underflow-resistant multiplication.
pub fn multiply_many_naive(distributions: &[&Array1<f64>]) -> ScaledArray {
  if distributions.is_empty() {
    return ScaledArray::empty(0);
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

  ScaledArray::new(result, log_scale)
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

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray {
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

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray {
    multiply_many_lazy_normalize(distributions)
  }
}

/// Log-scale aware multiplication with aggressive normalization.
pub struct AggressiveMultiply;

impl MultiplyAlgo for AggressiveMultiply {
  fn name(&self) -> &'static str {
    "aggressive-multiplication"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray {
    multiply_many(distributions)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_relative_eq;
  use ndarray::array;
  use std::iter::repeat_n;

  #[test]
  fn test_multiply_many_empty() {
    let result = multiply_many_lazy_normalize(&[]);
    assert_eq!(result.normalized.len(), 0);
    assert!(result.log_scale.is_infinite() && result.log_scale.is_sign_negative());
  }

  #[test]
  fn test_multiply_many_single() {
    let a = array![1.0, 4.0, 2.0];
    let result = multiply_many_lazy_normalize(&[&a]);

    assert_relative_eq!(result.normalized[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(result.normalized[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(result.normalized[2], 0.5, epsilon = 1e-10);
    assert_relative_eq!(result.log_scale, 4.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_multiply_many_two_identical() {
    let a = array![1.0, 2.0, 1.0];
    let result = multiply_many_lazy_normalize(&[&a, &a]);

    assert_relative_eq!(result.normalized[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(result.normalized[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(result.normalized[2], 0.25, epsilon = 1e-10);
    assert!(result.log_scale.is_finite());
  }

  #[test]
  fn test_multiply_many_underflow_resistance() {
    let n = 100;
    let small_val = 0.0001;
    let a = array![small_val];

    let refs: Vec<&Array1<f64>> = repeat_n(&a, n).collect();
    let result = multiply_many_lazy_normalize(&refs);

    assert!(!result.normalized.is_empty());
    assert!(result.log_scale.is_finite());

    let expected_log_scale = (n as f64) * small_val.ln();
    assert_relative_eq!(result.log_scale, expected_log_scale, epsilon = 1e-10);
  }

  #[test]
  fn test_multiply_many_handles_zero_max() {
    let a = array![0.0, 0.0, 0.0];
    let result = multiply_many_lazy_normalize(&[&a]);

    assert!(result.log_scale.is_infinite() && result.log_scale.is_sign_negative());
    assert_eq!(result.normalized, a);
  }

  #[test]
  fn test_lazy_vs_aggressive_same_result_for_few_factors() {
    let a = array![1.0, 2.0, 1.0];
    let b = array![0.5, 1.0, 0.5];

    let lazy = multiply_many_lazy_normalize(&[&a, &b]);
    let agg = multiply_many(&[&a, &b]);

    for i in 0..lazy.normalized.len() {
      assert_relative_eq!(lazy.normalized[i], agg.normalized[i], epsilon = 1e-10);
    }
    assert_relative_eq!(lazy.log_scale, agg.log_scale, epsilon = 1e-10);
  }

  #[test]
  fn test_lazy_normalize_fewer_operations() {
    let a = array![0.5, 1.0, 0.5];
    let b = array![0.5, 1.0, 0.5];
    let c = array![0.5, 1.0, 0.5];

    let result = multiply_many_lazy_normalize(&[&a, &b, &c]);

    assert_relative_eq!(result.normalized[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(result.normalized[0], 0.125, epsilon = 1e-10);
    assert_relative_eq!(result.normalized[2], 0.125, epsilon = 1e-10);
    assert_relative_eq!(result.log_scale, 0.0, epsilon = 1e-10);
  }
}
