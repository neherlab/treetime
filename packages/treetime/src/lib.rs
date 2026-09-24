pub mod alphabet;
pub mod ancestral;
pub mod cancel;
pub mod clock;
pub mod coalescent;
pub(crate) mod constants;
pub mod error;
pub mod gtr;
pub(crate) mod hacks;
pub mod homoplasy;
mod io;
pub mod mugration;
pub mod optimize;
pub mod partition;
pub mod progress;
pub mod prune;
pub mod reroot;
pub mod seq;
pub mod timetree;

#[cfg(test)]
pub(crate) mod test_utils;

#[cfg(test)]
mod graph;

pub use treetime_utils::{
  make_error, make_internal_error, make_internal_report, make_report, o, pretty_assert_abs_diff_eq,
  pretty_assert_neg_inf, pretty_assert_ulps_eq, vec_of_owned, vec_u8,
};

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
