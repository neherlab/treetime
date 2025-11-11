use crate::algos::algos::Algo;
use eyre::Report;
use ndarray::Array1;

/// Riemann convolution algorithm on uniform grids
pub fn convolve_riemann(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  let mut result = Array1::zeros(f_values.len() + g_values.len() - 1);

  for (i, &f_val) in f_values.iter().enumerate() {
    for (j, &g_val) in g_values.iter().enumerate() {
      result[i + j] += f_val * g_val * dx;
    }
  }

  Ok(result)
}

pub struct RiemannAlgo;

impl Algo for RiemannAlgo {
  fn name(&self) -> &'static str {
    "riemann"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_riemann(dx, f_values, g_values)
  }
}
