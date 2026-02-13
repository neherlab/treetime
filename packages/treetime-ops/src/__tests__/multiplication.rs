use crate::*;
use approx::assert_ulps_eq;
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

  assert_ulps_eq!(result.normalized[0], 0.25, max_ulps = 4);
  assert_ulps_eq!(result.normalized[1], 1.0, max_ulps = 4);
  assert_ulps_eq!(result.normalized[2], 0.5, max_ulps = 4);
  assert_ulps_eq!(result.log_scale, 4.0_f64.ln(), max_ulps = 4);
}

#[test]
fn test_multiply_many_two_identical() {
  let a = array![1.0, 2.0, 1.0];
  let result = multiply_many_lazy_normalize(&[&a, &a]);

  assert_ulps_eq!(result.normalized[0], 0.25, max_ulps = 4);
  assert_ulps_eq!(result.normalized[1], 1.0, max_ulps = 4);
  assert_ulps_eq!(result.normalized[2], 0.25, max_ulps = 4);
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
  assert_ulps_eq!(result.log_scale, expected_log_scale, max_ulps = 4);
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
    assert_ulps_eq!(lazy.normalized[i], agg.normalized[i], max_ulps = 4);
  }
  assert_ulps_eq!(lazy.log_scale, agg.log_scale, max_ulps = 4);
}

#[test]
fn test_lazy_normalize_fewer_operations() {
  let a = array![0.5, 1.0, 0.5];
  let b = array![0.5, 1.0, 0.5];
  let c = array![0.5, 1.0, 0.5];

  let result = multiply_many_lazy_normalize(&[&a, &b, &c]);

  assert_ulps_eq!(result.normalized[1], 1.0, max_ulps = 4);
  assert_ulps_eq!(result.normalized[0], 0.125, max_ulps = 4);
  assert_ulps_eq!(result.normalized[2], 0.125, max_ulps = 4);
  assert_ulps_eq!(result.log_scale, 0.0, max_ulps = 4);
}
