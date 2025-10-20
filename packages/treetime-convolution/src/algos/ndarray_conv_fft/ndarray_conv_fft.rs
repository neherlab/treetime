use crate::algos::algo_trait::Algo;
use crate::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

pub fn convolve_ndarray_conv_fft(f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  todo!()
}

/// Ndarray FFT-based convolution algorithm
pub struct NdarrayConvFftAlgo;

impl Algo for NdarrayConvFftAlgo {
  fn name(&self) -> &'static str {
    "ndarray"
  }

  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
    convolve_ndarray_conv_fft(f, g, x_grid)
  }
}
