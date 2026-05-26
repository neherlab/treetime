use crate::progress::ProgressSink;
use crate::{
  AncestralResult, ClockResult, MugrationResult, OptimizeResult, PruneResult, TimetreeResult, TreetimeAncestralArgs,
  TreetimeClockArgs, TreetimeMugrationArgs, TreetimeOptimizeArgs, TreetimePruneArgs, TreetimeTimetreeArgs,
};
use eyre::Report;
use treetime::commands::ancestral::run::run_ancestral_reconstruction;
use treetime::commands::clock::run::run_clock;
use treetime::commands::mugration::run::run_mugration;
use treetime::commands::optimize::run::run_optimize;
use treetime::commands::prune::run::run_prune;
use treetime::commands::timetree::run::run_timetree_estimation;

pub fn ancestral(args: &TreetimeAncestralArgs, progress: &dyn ProgressSink) -> Result<AncestralResult, Report> {
  run_ancestral_reconstruction(args, progress)
}

pub fn clock(args: &TreetimeClockArgs, progress: &dyn ProgressSink) -> Result<ClockResult, Report> {
  run_clock(args, progress)
}

pub fn timetree(args: &TreetimeTimetreeArgs, progress: &dyn ProgressSink) -> Result<TimetreeResult, Report> {
  run_timetree_estimation(args, progress)
}

pub fn mugration(args: &TreetimeMugrationArgs, progress: &dyn ProgressSink) -> Result<MugrationResult, Report> {
  run_mugration(args, progress)
}

pub fn optimize(args: &TreetimeOptimizeArgs, progress: &dyn ProgressSink) -> Result<OptimizeResult, Report> {
  run_optimize(args, progress)
}

pub fn prune(args: &TreetimePruneArgs, progress: &dyn ProgressSink) -> Result<PruneResult, Report> {
  run_prune(args, progress)
}
