pub mod algos;
pub mod grid_fn;
pub mod testing;

pub use algos::ndarray_conv::ndarray_conv::convolve_ndarray_conv as convolve;
pub use grid_fn::GridFn;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
