use crate::distribution::reference::convolution::{ndarray_convolve, riemann_convolve};
use crate::distribution::reference::convolution_test::traits::ConvAlgo;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

/// Riemann sum convolution algorithm
pub struct RiemannAlgo;

impl ConvAlgo for RiemannAlgo {
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
    riemann_convolve(f, g, x_grid)
  }

  fn name(&self) -> &'static str {
    "riemann"
  }
}

/// Ndarray FFT-based convolution algorithm
pub struct NdarrayAlgo;

impl ConvAlgo for NdarrayAlgo {
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
    ndarray_convolve(f, g, x_grid)
  }

  fn name(&self) -> &'static str {
    "ndarray"
  }
}
