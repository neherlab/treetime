use crate::distribution::reference::grid_fn::GridFn;
use crate::utils::container::minmax;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};

/// Convolution using ndarray-conv library (FFT-based) for uniform grids
pub fn ndarray_convolve(f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  assert!(
    is_uniform_grid(f.x()),
    "ndarray-conv requires uniform grid for function f"
  );
  assert!(
    is_uniform_grid(g.x()),
    "ndarray-conv requires uniform grid for function g"
  );
  assert!(is_uniform_grid(x_grid), "ndarray-conv requires uniform evaluation grid");

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

/// Convolution using Riemann sum integration for arbitrary grids
pub fn riemann_convolve(f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  let (f_min, f_max) = minmax(f.x().iter());
  let (g_min, g_max) = minmax(g.x().iter());

  let s_min = f_min - (g_max - g_min) * 0.1;
  let s_max = f_max + (g_max - g_min) * 0.1;
  let n_points = 2000;
  let ds = (s_max - s_min) / (n_points - 1) as f64;
  let s_grid = Array1::from_iter((0..n_points).map(|i| s_min + i as f64 * ds));

  let mut result = Array1::zeros(x_grid.len());
  for (i, &x) in x_grid.iter().enumerate() {
    let mut sum = 0.0;
    for &s in &s_grid {
      let f_val = f.interp(s).unwrap_or(0.0);
      let g_val = g.interp(x - s).unwrap_or(0.0);
      sum += f_val * g_val;
    }
    result[i] = sum * ds;
  }

  GridFn::new(x_grid.clone(), result)
}

/// Check if grid is uniform
fn is_uniform_grid(grid: &Array1<f64>) -> bool {
  if grid.len() < 2 {
    return true;
  }

  grid
    .iter()
    .tuple_windows()
    .map(|(a, b)| b - a)
    .tuple_windows()
    .all(|(d1, d2)| (d1 - d2).abs() <= 1e-12)
}
