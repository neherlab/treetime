pub mod algos;
pub mod analytical;
pub mod ops;
pub mod testing;

pub use treetime_grid::grid;
pub use treetime_grid::grid_fn;
pub use treetime_grid::grid_iter;
pub use treetime_grid::interp_nonuniform;
pub use treetime_grid::InterpElem;
pub use treetime_grid::{Grid, GridFn, GridFnF64};

pub use algos::ndarray_conv::ndarray_conv::convolve_ndarray_conv as convolve;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
