use eyre::Report;
use ndarray::Array1;
use treetime_ops::convolution::convolve_riemann;
use treetime_ops::ConvolveAlgo;

pub struct RiemannAlgo;

impl ConvolveAlgo for RiemannAlgo {
  fn name(&self) -> &'static str {
    "riemann"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_riemann(dx, f_values, g_values)
  }
}
