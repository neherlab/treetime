use crate::algos::algo_trait::Algo;
use crate::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use treetime_utils::ndarray::is_uniform_grid;

/// Convolution using ndarray_conv library for uniform grids
pub fn convolve_ndarray_conv(f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  debug_assert!(is_uniform_grid(f.x()), "grid must be uniform");
  debug_assert!(is_uniform_grid(g.x()), "grid must be uniform");
  debug_assert!(is_uniform_grid(x_grid), "grid must be uniform");

  let raw_conv = f.y().conv(g.y(), ConvMode::Full, PaddingMode::Zeros)?;
  let dx = f.x()[1] - f.x()[0];
  let conv_scaled = &raw_conv * dx;

  let result_len = raw_conv.len();
  let result_min = f.x()[0] + g.x()[0];
  let result_x = Array1::from_iter((0..result_len).map(|i| result_min + i as f64 * dx));

  if result_x.len() == x_grid.len() {
    GridFn::new(result_x, conv_scaled)
  } else {
    let temp_fn = GridFn::new(result_x, conv_scaled)?;
    let result_vals = temp_fn.interp_many(x_grid)?;
    GridFn::new(x_grid.clone(), result_vals)
  }
}

pub struct NdarrayAlgo;

impl Algo for NdarrayAlgo {
  fn name(&self) -> &'static str {
    "ndarray"
  }

  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
    convolve_ndarray_conv(f, g, x_grid)
  }
}
