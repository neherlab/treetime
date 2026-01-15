pub mod convolution;
pub mod multiplication;
pub mod traits;

pub use convolution::{FftConvolve, NdarrayConvolve, RiemannConvolve, convolve, convolve_fft, convolve_riemann};
pub use multiplication::{
  LogScaleMultiply, PointwiseMultiply, multiply_many, multiply_many_lazy_normalize, multiply_many_naive,
};
pub use traits::{ConvolveAlgo, MultiplyAlgo};

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
