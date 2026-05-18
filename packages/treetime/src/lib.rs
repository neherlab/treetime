pub mod alphabet;
pub mod ancestral;
pub mod clock;
pub mod coalescent;
#[cfg(feature = "clap")]
pub mod commands;
pub mod constants;
pub mod graph;
pub mod gtr;
pub mod hacks;
pub mod io;
pub mod optimize;
pub mod partition;
pub mod seq;
pub mod timetree;

#[cfg(test)]
pub mod test_utils;

pub use treetime_utils::{
  make_error, make_internal_error, make_internal_report, make_report, o, pretty_assert_abs_diff_eq,
  pretty_assert_neg_inf, pretty_assert_ulps_eq, vec_of_owned, vec_u8,
};

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
