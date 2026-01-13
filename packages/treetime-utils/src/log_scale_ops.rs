use ndarray::Array1;

fn array_max(arr: &Array1<f64>) -> f64 {
  arr.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Multiply multiple distributions in log-scale to avoid underflow.
///
/// Returns (normalized_shape, log_scale) where the full result is:
/// `normalized_shape * exp(log_scale)`
///
/// Invariant: max(normalized_shape) = 1.0
///
/// After each pairwise multiplication, the result is normalized (max = 1.0)
/// and the scale factor is accumulated in log-space. This prevents underflow
/// when multiplying many small probability values.
pub fn log_scale_multiply_many(distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
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

/// Multiply distributions naively without log-scale tracking.
///
/// Returns (normalized_shape, log_scale) where the full result is:
/// `normalized_shape * exp(log_scale)`
///
/// Warning: This can underflow to zero for many small distributions.
/// Use `log_scale_multiply_many` for underflow-resistant multiplication.
pub fn naive_multiply_many(distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
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

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_relative_eq;
  use ndarray::array;

  #[test]
  fn test_log_scale_multiply_empty() {
    let (shape, log_scale) = log_scale_multiply_many(&[]);
    assert_eq!(shape.len(), 0);
    assert!(log_scale.is_infinite() && log_scale.is_sign_negative());
  }

  #[test]
  fn test_log_scale_multiply_single() {
    let a = array![1.0, 4.0, 2.0];
    let (shape, log_scale) = log_scale_multiply_many(&[&a]);

    assert_relative_eq!(shape[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(shape[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(shape[2], 0.5, epsilon = 1e-10);
    assert_relative_eq!(log_scale, 4.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_log_scale_multiply_two_identical() {
    let a = array![1.0, 2.0, 1.0];
    let (shape, log_scale) = log_scale_multiply_many(&[&a, &a]);

    assert_relative_eq!(shape[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(shape[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(shape[2], 0.25, epsilon = 1e-10);
    assert!(log_scale.is_finite());
  }

  #[test]
  fn test_log_scale_multiply_underflow_resistance() {
    let n = 100;
    let small_val = 0.0001;
    let a = array![small_val];

    let refs: Vec<&Array1<f64>> = std::iter::repeat_n(&a, n).collect();
    let (shape, log_scale) = log_scale_multiply_many(&refs);

    assert!(!shape.is_empty());
    assert!(log_scale.is_finite());

    let expected_log_scale = (n as f64) * small_val.ln();
    assert_relative_eq!(log_scale, expected_log_scale, epsilon = 1e-10);
  }

  #[test]
  fn test_naive_multiply_empty() {
    let (shape, log_scale) = naive_multiply_many(&[]);
    assert_eq!(shape.len(), 0);
    assert!(log_scale.is_infinite() && log_scale.is_sign_negative());
  }

  #[test]
  fn test_naive_multiply_single() {
    let a = array![1.0, 4.0, 2.0];
    let (shape, log_scale) = naive_multiply_many(&[&a]);

    assert_relative_eq!(shape[0], 0.25, epsilon = 1e-10);
    assert_relative_eq!(shape[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(shape[2], 0.5, epsilon = 1e-10);
    assert_relative_eq!(log_scale, 4.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_log_scale_handles_zero_max() {
    let a = array![0.0, 0.0, 0.0];
    let (shape, log_scale) = log_scale_multiply_many(&[&a]);

    assert!(log_scale.is_infinite() && log_scale.is_sign_negative());
    assert_eq!(shape, a);
  }
}
