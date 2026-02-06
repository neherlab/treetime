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
