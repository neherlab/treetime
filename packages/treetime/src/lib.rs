pub mod alphabet;
pub mod cli;
pub mod commands;
pub mod constants;
pub mod distribution;
pub mod graph;
pub mod gtr;
pub mod hacks;
pub mod io;
pub mod representation;
pub mod seq;

pub use treetime_graph;
pub use treetime_utils::{
  make_error, make_internal_error, make_internal_report, make_report, o, pretty_assert_abs_diff_eq,
  pretty_assert_ulps_eq, vec_of_owned, vec_u8,
};

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::global_init::global_init;

  #[ctor]
  fn init() {
    global_init();
  }
}
