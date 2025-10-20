use eyre::Report;
use ndarray::Array1;

pub trait Algo: Send + Sync {
  fn name(&self) -> &'static str;

  fn convolve(
    &self,
    input_grid: &Array1<f64>,
    f_values: &Array1<f64>,
    g_values: &Array1<f64>,
    output_grid: &Array1<f64>,
  ) -> Result<Array1<f64>, Report>;
}
