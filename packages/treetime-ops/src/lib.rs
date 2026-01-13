pub mod convolution;
pub mod multiplication;
pub mod traits;

pub use convolution::{convolve, convolve_fft, convolve_riemann};
pub use multiplication::{multiply_many, multiply_many_lazy_normalize};
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
