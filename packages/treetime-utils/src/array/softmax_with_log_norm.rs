use crate::cast::ToF64Exact;
use ndarray::ArrayViewMut1;
use ndarray_stats::QuantileExt;
use ndarray_stats::errors::MinMaxError;

pub fn softmax_with_log_norm_inplace(mut log_vec: ArrayViewMut1<'_, f64>) -> f64 {
  let n = log_vec.len();

  let max_val = match log_vec.max() {
    Ok(&max) => max,
    Err(MinMaxError::EmptyInput) => return f64::NEG_INFINITY,
    Err(MinMaxError::UndefinedOrder) => {
      log_vec.fill(f64::NAN);
      return f64::NAN;
    },
  };

  if max_val == f64::INFINITY {
    let inf_count = log_vec.iter().filter(|&&v| v == f64::INFINITY).count();
    let prob = 1.0 / inf_count.to_f64_exact();
    log_vec.mapv_inplace(|v| if v == f64::INFINITY { prob } else { 0.0 });
    return f64::INFINITY;
  }

  if !max_val.is_finite() {
    log_vec.fill(1.0 / n.to_f64_exact());
    return f64::NEG_INFINITY;
  }

  log_vec.mapv_inplace(|v| (v - max_val).exp());
  let shifted_sum = log_vec.sum();

  let log_norm = max_val + shifted_sum.ln();

  log_vec /= shifted_sum;

  log_norm
}
