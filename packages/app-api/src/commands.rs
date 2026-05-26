use crate::progress::ProgressSink;
use eyre::Report;
use treetime::commands::ancestral::args::TreetimeAncestralArgs;
use treetime::commands::ancestral::run::run_ancestral_reconstruction;
use treetime::commands::clock::args::TreetimeClockArgs;
use treetime::commands::clock::run::run_clock;
use treetime::commands::mugration::args::TreetimeMugrationArgs;
use treetime::commands::mugration::run::run_mugration;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::optimize::run::run_optimize;
use treetime::commands::prune::args::TreetimePruneArgs;
use treetime::commands::prune::run::run_prune;
use treetime::commands::timetree::args::TreetimeTimetreeArgs;
use treetime::commands::timetree::run::run_timetree_estimation;

pub fn ancestral(args: &TreetimeAncestralArgs, _progress: &dyn ProgressSink) -> Result<(), Report> {
  drop(run_ancestral_reconstruction(args)?);
  Ok(())
}

pub fn clock(args: &TreetimeClockArgs, _progress: &dyn ProgressSink) -> Result<(), Report> {
  drop(run_clock(args)?);
  Ok(())
}

pub fn timetree(args: &TreetimeTimetreeArgs, _progress: &dyn ProgressSink) -> Result<(), Report> {
  drop(run_timetree_estimation(args)?);
  Ok(())
}

pub fn mugration(args: &TreetimeMugrationArgs, _progress: &dyn ProgressSink) -> Result<(), Report> {
  drop(run_mugration(args)?);
  Ok(())
}

pub fn optimize(args: &TreetimeOptimizeArgs, _progress: &dyn ProgressSink) -> Result<(), Report> {
  drop(run_optimize(args)?);
  Ok(())
}

pub fn prune(args: &TreetimePruneArgs, _progress: &dyn ProgressSink) -> Result<(), Report> {
  drop(run_prune(args)?);
  Ok(())
}
