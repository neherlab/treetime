pub mod commands;
pub mod datasets;
pub mod pipelines;
pub mod progress;

pub use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
pub use crate::commands::ancestral::result::AncestralResult;
pub use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
pub use crate::commands::clock::run::ClockResult;
pub use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
pub use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
pub use crate::commands::optimize::result::OptimizeResult;
pub use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
pub use crate::commands::prune::result::PruneResult;
pub use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
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
