pub mod commands;
pub mod progress;

pub use treetime::commands::ancestral::args::TreetimeAncestralArgs;
pub use treetime::commands::clock::args::TreetimeClockArgs;
pub use treetime::commands::mugration::args::TreetimeMugrationArgs;
pub use treetime::commands::optimize::args::TreetimeOptimizeArgs;
pub use treetime::commands::prune::args::TreetimePruneArgs;
pub use treetime::commands::timetree::args::TreetimeTimetreeArgs;
pub use treetime::version;

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
