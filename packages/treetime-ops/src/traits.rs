use crate::ScaledArray;
use eyre::Report;
use ndarray::Array1;

pub trait ConvolveAlgo: Send + Sync {
  fn name(&self) -> &'static str;

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report>;
}

pub trait MultiplyAlgo: Send + Sync {
  fn name(&self) -> &'static str;

  fn multiply(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Array1<f64>;

  fn multiply_many(&self, distributions: &[&Array1<f64>]) -> ScaledArray;
}
