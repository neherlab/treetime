use eyre::Report;
use ndarray::Array1;
use treetime_ops::convolution::convolve_fft;
use treetime_ops::ConvolveAlgo;

pub struct NdarrayConvFftAlgo;

impl ConvolveAlgo for NdarrayConvFftAlgo {
  fn name(&self) -> &'static str {
    "ndarray-conv-fft"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_fft(dx, f_values, g_values)
  }
}
