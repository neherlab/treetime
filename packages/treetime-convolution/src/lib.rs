pub mod algos;
pub mod grid_fn;
pub mod testing;

use num_traits::{Num, NumCast};
use std::fmt::Debug;

pub trait InterpElem: Num + NumCast + Debug + Send + PartialOrd + Copy {}

impl InterpElem for f64 {}

pub use algos::ndarray_conv::ndarray_conv::convolve_ndarray_conv as convolve;
pub use grid_fn::GridFn;
pub type GridFnF64 = GridFn<f64>;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
