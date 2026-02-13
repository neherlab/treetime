use ndarray::Array1;
use serde::{Deserialize, Serialize};

pub mod convolution;
pub mod multiplication;
pub mod traits;

#[cfg(test)]
mod __tests__;

pub use convolution::{FftConvolve, NdarrayConvolve, RiemannConvolve, convolve, convolve_fft, convolve_riemann};
pub use multiplication::{
  AggressiveMultiply, LogScaleMultiply, PointwiseMultiply, multiply_many, multiply_many_lazy_normalize,
  multiply_many_naive,
};
pub use traits::{ConvolveAlgo, MultiplyAlgo};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScaledArray {
  pub normalized: Array1<f64>,
  pub log_scale: f64,
}

impl ScaledArray {
  pub fn new(normalized: Array1<f64>, log_scale: f64) -> Self {
    Self { normalized, log_scale }
  }

  pub fn empty(len: usize) -> Self {
    Self {
      normalized: Array1::zeros(len),
      log_scale: f64::NEG_INFINITY,
    }
  }
}

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
