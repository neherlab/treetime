use crate::ScaledArray;
use eyre::Report;
use ndarray::Array1;

/// Trait for convolution algorithms operating on discrete probability distributions.
pub trait ConvolveAlgo: Send + Sync {
  fn name(&self) -> &'static str;

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report>;
}

/// Trait for multiplication algorithms operating on discrete probability distributions.
///
/// Unlike convolution, point-wise multiplication does not depend on grid spacing.
pub trait MultiplyAlgo: Send + Sync {
  fn name(&self) -> &'static str;

  /// Multiply two distributions element-wise.
  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64>;

  /// Multiply N distributions, returning normalized shape and log-scale.
  ///
  /// Returns `ScaledArray { normalized, log_scale }` where the full result is
  /// `normalized * exp(log_scale)`.
  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray;
}
