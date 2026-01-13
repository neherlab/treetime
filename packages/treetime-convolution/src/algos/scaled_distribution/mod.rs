use ndarray::Array1;
use treetime_ops::multiplication::{multiply_many_lazy_normalize, multiply_many_naive};
use treetime_ops::MultiplyAlgo;

/// Naive point-wise multiplication algorithm.
///
/// Multiplies distributions element-wise without any log-scale handling.
/// This is the simplest possible implementation for testing purposes.
pub struct NaiveMultiplicationAlgo;

impl MultiplyAlgo for NaiveMultiplicationAlgo {
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

/// Log-scale aware multiplication algorithm with lazy normalization.
///
/// Multiplies distributions while tracking scale in log-space to avoid underflow.
/// Uses lazy normalization - only normalizes when approaching underflow threshold,
/// improving both performance and accuracy compared to aggressive normalization.
pub struct LogScaleMultiplicationAlgo;

impl MultiplyAlgo for LogScaleMultiplicationAlgo {
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
