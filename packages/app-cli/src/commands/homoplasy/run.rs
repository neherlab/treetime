use crate::commands::homoplasy::args::TreetimeHomoplasyArgs;
use crate::commands::homoplasy::result::HomoplasyResult;
use eyre::Report;
use treetime::cancel::Cancel;
use treetime::homoplasy::pipeline::{self, HomoplasyInput, HomoplasyParams};
use treetime::progress::ProgressSink;

pub(crate) fn run_homoplasy(
  _: &TreetimeHomoplasyArgs,
  cancel: &dyn Cancel,
  _: &dyn ProgressSink,
) -> Result<HomoplasyResult, Report> {
  cancel.check()?;
  pipeline::run(&HomoplasyParams, HomoplasyInput).map_err(|err| err.into_report())?;
  Ok(HomoplasyResult)
}
