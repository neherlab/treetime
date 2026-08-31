pub mod commands;
pub mod datasets;
pub mod pipelines;
pub mod progress;

pub use treetime::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
pub use treetime::commands::ancestral::result::AncestralResult;
pub use treetime::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
pub use treetime::commands::clock::run::ClockResult;
pub use treetime::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
pub use treetime::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
pub use treetime::commands::optimize::result::OptimizeResult;
pub use treetime::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
pub use treetime::commands::prune::result::PruneResult;
pub use treetime::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
pub use treetime::commands::timetree::result::TimetreeResult;
pub use treetime::mugration::result::MugrationResult;

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
