use crate::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

pub trait Algo: Send + Sync {
  fn name(&self) -> &'static str;

  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;
}
