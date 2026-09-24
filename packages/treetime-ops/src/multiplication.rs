use crate::ScaledArray;
use crate::traits::MultiplyAlgo;
use ndarray::Array1;
use treetime_utils::array::ndarray::max_or;

const UNDERFLOW_THRESHOLD: f64 = 1e-100;

pub struct PointwiseMultiply;

impl MultiplyAlgo for PointwiseMultiply {
  fn name(&self) -> &'static str {
    "pointwise"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray {
    multiply_many_naive(distributions)
  }
}

fn multiply_many_naive(distributions: &[&Array1<f64>]) -> ScaledArray {
  if distributions.is_empty() {
    return ScaledArray::empty(0);
  }

  let n = distributions[0].len();
  let mut result = Array1::ones(n);

  for dist in distributions {
    result = &result * *dist;
  }

  let max_val = max_or(&result, f64::NEG_INFINITY);
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

pub struct LogScaleMultiply;

impl MultiplyAlgo for LogScaleMultiply {
  fn name(&self) -> &'static str {
    "log-scale"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray {
    multiply_many_lazy_normalize(distributions)
  }
}

pub(crate) fn multiply_many_lazy_normalize(distributions: &[&Array1<f64>]) -> ScaledArray {
  if distributions.is_empty() {
    return ScaledArray::empty(0);
  }

  if distributions.len() == 1 {
    return normalize_and_extract_scale(distributions[0]);
  }

  let n = distributions[0].len();
  let mut product = distributions[0].clone();
  let mut accumulated_log_scale = 0.0;

  let max_val = max_or(&product, f64::NEG_INFINITY);
  if max_val <= 0.0 || !max_val.is_finite() {
    return ScaledArray::empty(n);
  }

  for dist in distributions.iter().skip(1) {
    product = &product * *dist;

    let max_val = max_or(&product, f64::NEG_INFINITY);
    if max_val <= 0.0 || !max_val.is_finite() {
      return ScaledArray::empty(n);
    }

    if max_val < UNDERFLOW_THRESHOLD {
      accumulated_log_scale += max_val.ln();
      product.mapv_inplace(|v| v / max_val);
    }
  }

  let max_val = max_or(&product, f64::NEG_INFINITY);
  if max_val <= 0.0 || !max_val.is_finite() {
    return ScaledArray::empty(n);
  }
  accumulated_log_scale += max_val.ln();
  product.mapv_inplace(|v| v / max_val);

  ScaledArray::new(product, accumulated_log_scale)
}

pub struct AggressiveMultiply;

impl MultiplyAlgo for AggressiveMultiply {
  fn name(&self) -> &'static str {
    "aggressive"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray {
    multiply_many(distributions)
  }
}

pub(crate) fn multiply_many(distributions: &[&Array1<f64>]) -> ScaledArray {
  if distributions.is_empty() {
    return ScaledArray::empty(0);
  }

  if distributions.len() == 1 {
    return normalize_and_extract_scale(distributions[0]);
  }

  let n = distributions[0].len();
  let mut accumulated_log_scale = 0.0;
  let mut normalized_result = distributions[0].clone();

  let max_val = max_or(&normalized_result, f64::NEG_INFINITY);
  if max_val <= 0.0 || !max_val.is_finite() {
    return ScaledArray::empty(n);
  }
  accumulated_log_scale += max_val.ln();
  normalized_result.mapv_inplace(|v| v / max_val);

  for dist in distributions.iter().skip(1) {
    normalized_result = &normalized_result * *dist;

    let max_val = max_or(&normalized_result, f64::NEG_INFINITY);
    if max_val <= 0.0 || !max_val.is_finite() {
      return ScaledArray::empty(n);
    }

    accumulated_log_scale += max_val.ln();
    normalized_result.mapv_inplace(|v| v / max_val);
  }

  ScaledArray::new(normalized_result, accumulated_log_scale)
}

fn normalize_and_extract_scale(arr: &Array1<f64>) -> ScaledArray {
  let max_val = max_or(arr, f64::NEG_INFINITY);
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
