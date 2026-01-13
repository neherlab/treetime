use eyre::Report;
use ndarray::Array1;
use treetime_ops::convolution::convolve;
use treetime_ops::ConvolveAlgo;

pub use treetime_ops::convolution::convolve as convolve_ndarray_conv;

pub struct NdarrayAlgo;

impl ConvolveAlgo for NdarrayAlgo {
  fn name(&self) -> &'static str {
    "ndarray-conv"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve(dx, f_values, g_values)
  }
}
