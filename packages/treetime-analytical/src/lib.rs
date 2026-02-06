pub mod exponential;
pub mod gaussian;
pub mod gaussian_exponential;
pub mod validation;

pub use exponential::{exponential_convolution, exponential_convolution_grid, exponential_pdf, exponential_pdf_grid};
pub use gaussian::{
  GaussianParams, gaussian_convolution, gaussian_convolution_pdf, gaussian_convolution_pdf_grid, gaussian_evaluate,
  gaussian_pdf, gaussian_pdf_grid, gaussian_product, gaussian_product_params,
};
pub use gaussian_exponential::{gaussian_exponential_convolution, gaussian_exponential_convolution_grid};

#[cfg(test)]
mod gaussian_tests;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
