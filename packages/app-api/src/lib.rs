pub mod commands;
pub mod progress;
pub mod version;

pub use treetime::commands::ancestral::args::TreetimeAncestralArgs;
pub use treetime::commands::ancestral::result::AncestralResult;
pub use treetime::commands::clock::args::TreetimeClockArgs;
pub use treetime::commands::clock::run::ClockResult;
pub use treetime::commands::mugration::args::TreetimeMugrationArgs;
pub use treetime::commands::optimize::args::TreetimeOptimizeArgs;
pub use treetime::commands::optimize::result::OptimizeResult;
pub use treetime::commands::prune::args::TreetimePruneArgs;
pub use treetime::commands::prune::result::PruneResult;
pub use treetime::commands::timetree::args::TreetimeTimetreeArgs;
<<<<<<< HEAD
pub use treetime::version;
=======
pub use treetime::commands::timetree::result::TimetreeResult;
pub use treetime::mugration::result::MugrationResult;
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)

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
