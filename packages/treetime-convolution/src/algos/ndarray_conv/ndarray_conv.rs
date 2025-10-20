use crate::algos::algos::Algo;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use ndarray_interp::interp1d::Interp1DBuilder;
use treetime_utils::ndarray::is_uniform_grid;

/// Convolution using ndarray_conv library for uniform grids
pub fn convolve_ndarray_conv(
  input_grid: &Array1<f64>,
  f_values: &Array1<f64>,
  g_values: &Array1<f64>,
  output_grid: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  debug_assert!(is_uniform_grid(input_grid), "input_grid must be uniform");
  debug_assert!(is_uniform_grid(output_grid), "output_grid must be uniform");

  let grid_spacing = input_grid[1] - input_grid[0];

  let discrete_conv = f_values.conv(g_values, ConvMode::Full, PaddingMode::Zeros)?;
  let continuous_conv = &discrete_conv * grid_spacing;

  let conv_result_len = discrete_conv.len();
  let conv_result_min = input_grid[0] + input_grid[0];
  let conv_result_max = input_grid[input_grid.len() - 1] + input_grid[input_grid.len() - 1];
  let conv_result_grid = Array1::linspace(conv_result_min, conv_result_max, conv_result_len);

  let temp_interp = Interp1DBuilder::new(continuous_conv).x(conv_result_grid).build()?;

  let mut result = Array1::zeros(output_grid.len());
  for (i, &x) in output_grid.iter().enumerate() {
    result[i] = temp_interp.interp_scalar(x).unwrap_or(0.0);
  }

  Ok(result)
}

pub struct NdarrayAlgo;

impl Algo for NdarrayAlgo {
  fn name(&self) -> &'static str {
    "ndarray-conv"
  }

  fn convolve(
    &self,
    input_grid: &Array1<f64>,
    f_values: &Array1<f64>,
    g_values: &Array1<f64>,
    output_grid: &Array1<f64>,
  ) -> Result<Array1<f64>, Report> {
    convolve_ndarray_conv(input_grid, f_values, g_values, output_grid)
  }
}
