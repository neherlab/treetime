pub mod grid;
pub mod grid_fn;
pub mod grid_iter;
pub mod interp_nonuniform;

use num_traits::{Num, NumCast};
use std::fmt::Debug;

pub trait InterpElem: Num + NumCast + Debug + Send + PartialOrd + Copy {}

impl InterpElem for f64 {}

pub use grid::Grid;
pub use grid_fn::GridFn;
pub type GridFnF64 = GridFn<f64>;

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
