pub mod gaussian;

pub use gaussian::{GaussianParams, gaussian_product, gaussian_product_params};

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
