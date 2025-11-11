use crate::algos::algos::Algo;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};

/// Convolution using ndarray_conv library for uniform grids
/// takes as input the grid spacing, and arrays with function values
pub fn convolve_ndarray_conv(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  let discrete_conv = f_values.conv(g_values, ConvMode::Full, PaddingMode::Zeros)?;
  let continuous_conv = &discrete_conv * dx;
  Ok(continuous_conv)
}

pub struct NdarrayAlgo;

impl Algo for NdarrayAlgo {
  fn name(&self) -> &'static str {
    "ndarray-conv"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_ndarray_conv(dx, f_values, g_values)
  }
}
