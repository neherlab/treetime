use ndarray::{Array1, ArrayView1};
use ndarray_stats::QuantileExt;
use ndarray_stats::errors::MinMaxError;

#[allow(
  clippy::as_conversions,
  reason = "element counts to f64 for uniform normalization; exact for any realistic length"
)]
pub fn softmax_with_log_norm(log_vec: ArrayView1<'_, f64>) -> (Array1<f64>, f64) {
  let n = log_vec.len();

  let max_val = match log_vec.max() {
    Ok(&max) => max,
    Err(MinMaxError::EmptyInput) => return (Array1::zeros(0), f64::NEG_INFINITY),
    Err(MinMaxError::UndefinedOrder) => return (Array1::from_elem(n, f64::NAN), f64::NAN),
  };

  if max_val == f64::INFINITY {
    let inf_count = log_vec.iter().filter(|&&v| v == f64::INFINITY).count();
    let prob = 1.0 / inf_count as f64;
    return (
      log_vec.mapv(|v| if v == f64::INFINITY { prob } else { 0.0 }),
      f64::INFINITY,
    );
  }

  if !max_val.is_finite() {
    return (Array1::from_elem(n, 1.0 / n as f64), f64::NEG_INFINITY);
  }

  let shifted = log_vec.mapv(|v| (v - max_val).exp());
  let shifted_sum = shifted.sum();

  let log_norm = max_val + shifted_sum.ln();

  let normalized = shifted / shifted_sum;

  (normalized, log_norm)
}
