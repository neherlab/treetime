use crate::algos::algos::MultiplyAlgo;
use eyre::Report;
use ndarray::Array1;

/// Naive point-wise multiplication algorithm.
///
/// Multiplies distributions element-wise without any log-scale handling.
/// This is the simplest possible implementation for testing purposes.
pub struct NaiveMultiplicationAlgo;

impl MultiplyAlgo for NaiveMultiplicationAlgo {
  fn name(&self) -> &'static str {
    "naive-multiplication"
  }

  fn multiply(&self, _dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(f_values * g_values)
  }

  fn multiply_many(&self, _dx: f64, distributions: &[&Array1<f64>]) -> Result<(Array1<f64>, f64), Report> {
    if distributions.is_empty() {
      return Ok((Array1::zeros(0), f64::NEG_INFINITY));
    }

    let n = distributions[0].len();
    let mut result = Array1::ones(n);

    for dist in distributions {
      result = &result * *dist;
    }

    let max_val = result.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let log_scale = if max_val > 0.0 && max_val.is_finite() {
      max_val.ln()
    } else {
      f64::NEG_INFINITY
    };

    if max_val > 0.0 && max_val.is_finite() {
      result.mapv_inplace(|v| v / max_val);
    }

    Ok((result, log_scale))
  }
}

/// Log-scale aware multiplication algorithm.
///
/// Multiplies distributions while tracking scale in log-space to avoid underflow.
/// This mimics the behavior of ScaledDistribution multiplication.
pub struct LogScaleMultiplicationAlgo;

impl MultiplyAlgo for LogScaleMultiplicationAlgo {
  fn name(&self) -> &'static str {
    "log-scale-multiplication"
  }

  fn multiply(&self, _dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let product = f_values * g_values;
    let max_val = product.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    if max_val <= 0.0 || !max_val.is_finite() {
      return Ok(Array1::zeros(f_values.len()));
    }

    Ok(product)
  }

  fn multiply_many(&self, _dx: f64, distributions: &[&Array1<f64>]) -> Result<(Array1<f64>, f64), Report> {
    if distributions.is_empty() {
      return Ok((Array1::zeros(0), f64::NEG_INFINITY));
    }

    if distributions.len() == 1 {
      let max_val = distributions[0].iter().copied().fold(f64::NEG_INFINITY, f64::max);
      let log_scale = if max_val > 0.0 && max_val.is_finite() {
        max_val.ln()
      } else {
        f64::NEG_INFINITY
      };
      let normalized = if max_val > 0.0 && max_val.is_finite() {
        distributions[0].mapv(|v| v / max_val)
      } else {
        distributions[0].clone()
      };
      return Ok((normalized, log_scale));
    }

    let n = distributions[0].len();
    let mut accumulated_log_scale = 0.0;
    let mut normalized_result = distributions[0].clone();

    let max_val = normalized_result.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_val <= 0.0 || !max_val.is_finite() {
      return Ok((Array1::zeros(n), f64::NEG_INFINITY));
    }
    accumulated_log_scale += max_val.ln();
    normalized_result.mapv_inplace(|v| v / max_val);

    for dist in distributions.iter().skip(1) {
      normalized_result = &normalized_result * *dist;

      let max_val = normalized_result.iter().copied().fold(f64::NEG_INFINITY, f64::max);
      if max_val <= 0.0 || !max_val.is_finite() {
        return Ok((Array1::zeros(n), f64::NEG_INFINITY));
      }

      accumulated_log_scale += max_val.ln();
      normalized_result.mapv_inplace(|v| v / max_val);
    }

    Ok((normalized_result, accumulated_log_scale))
  }
}
