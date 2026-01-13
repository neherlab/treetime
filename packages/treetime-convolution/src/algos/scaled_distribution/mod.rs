use crate::algos::algos::MultiplyAlgo;
use ndarray::Array1;
use treetime_utils::log_scale_ops::{log_scale_multiply_many, naive_multiply_many};

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
    naive_multiply_many(distributions)
  }
}

/// Log-scale aware multiplication algorithm.
///
/// Multiplies distributions while tracking scale in log-space to avoid underflow.
/// After each pairwise multiplication, the result is normalized (max = 1.0)
/// and the scale factor is accumulated in log-space.
pub struct LogScaleMultiplicationAlgo;

impl MultiplyAlgo for LogScaleMultiplicationAlgo {
  fn name(&self) -> &'static str {
    "log-scale-multiplication"
  }

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64> {
    f_values * g_values
  }

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> (Array1<f64>, f64) {
    log_scale_multiply_many(distributions)
  }
}
