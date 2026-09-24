pub mod exponential;
pub mod gaussian;
pub mod gaussian_exponential;
pub mod validation;

pub use exponential::{exponential_convolution_grid, exponential_pdf_grid};
pub use gaussian::{GaussianParams, gaussian_convolution_pdf_grid, gaussian_pdf_grid, gaussian_product};
pub use gaussian_exponential::gaussian_exponential_convolution_grid;

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
