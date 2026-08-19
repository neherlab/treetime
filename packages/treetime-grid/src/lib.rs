pub mod boundary_behavior;
pub mod grid;
pub mod grid_edge;
pub mod grid_fn;
pub mod grid_iter;
pub mod hard_approach_law;
pub mod interp_nonuniform;
pub mod piecewise_constant_fn;
pub mod piecewise_fn;
pub mod piecewise_linear_fn;
pub mod soft_tail_law;

use num_traits::{Num, NumCast};
use std::fmt::Debug;

pub trait InterpElem: Num + NumCast + Debug + Send + PartialOrd + Copy {}

impl InterpElem for f64 {}

pub use boundary_behavior::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS};
pub use grid::Grid;
pub use grid_edge::GridEdge;
pub use grid_fn::GridFn;
pub use hard_approach_law::{Approach, HardApproachLaw, Side};
pub use soft_tail_law::SoftTailLaw;
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
