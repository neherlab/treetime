use crate::algos::algo_trait::Algo;
use eyre::Report;
use ndarray::Array1;

pub fn convolve_ndarray_conv_fft(
  _input_grid: &Array1<f64>,
  _f_values: &Array1<f64>,
  _g_values: &Array1<f64>,
  _output_grid: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  todo!()
}

/// Ndarray FFT-based convolution algorithm
pub struct NdarrayConvFftAlgo;

impl Algo for NdarrayConvFftAlgo {
  fn name(&self) -> &'static str {
    "ndarray_fft"
  }

  fn convolve(
    &self,
    input_grid: &Array1<f64>,
    f_values: &Array1<f64>,
    g_values: &Array1<f64>,
    output_grid: &Array1<f64>,
  ) -> Result<Array1<f64>, Report> {
    convolve_ndarray_conv_fft(input_grid, f_values, g_values, output_grid)
  }
}
