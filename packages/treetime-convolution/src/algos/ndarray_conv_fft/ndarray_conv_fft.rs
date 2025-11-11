use crate::algos::algos::Algo;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvFFTExt, ConvMode, PaddingMode};

/// Convolution using ndarray_conv library with FFT for uniform grids.
/// Takes as input the grid spacing, and arrays with function values
pub fn convolve_ndarray_conv_fft(
  dx: f64,
  f_values: &Array1<f64>,
  g_values: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  let discrete_conv = f_values.conv_fft(g_values, ConvMode::Full, PaddingMode::Zeros)?;
  let continuous_conv = &discrete_conv * dx;
  Ok(continuous_conv)
}

/// Ndarray FFT-based convolution algorithm
pub struct NdarrayConvFftAlgo;

impl Algo for NdarrayConvFftAlgo {
  fn name(&self) -> &'static str {
    "ndarray-conv-fft"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_ndarray_conv_fft(dx, f_values, g_values)
  }
}
