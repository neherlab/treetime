use crate::algos::algo_trait::Algo;
use crate::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;
use treetime_utils::container::minmax;

/// Convolution using Riemann sum integration for arbitrary grids
pub fn convolve_riemann(f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
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

pub struct RiemannAlgo;

impl Algo for RiemannAlgo {
  fn name(&self) -> &'static str {
    "riemann"
  }

  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
    convolve_riemann(f, g, x_grid)
  }
}
