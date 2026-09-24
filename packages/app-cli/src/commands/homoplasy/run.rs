use crate::commands::homoplasy::args::TreetimeHomoplasyArgs;
use crate::commands::homoplasy::result::HomoplasyResult;
use eyre::Report;
use treetime::homoplasy::pipeline::{self, HomoplasyInput, HomoplasyParams};

pub(crate) fn run_homoplasy(_: TreetimeHomoplasyArgs) -> Result<HomoplasyResult, Report> {
  pipeline::run(&HomoplasyParams, HomoplasyInput).map_err(|err| err.into_report())?;
  Ok(HomoplasyResult)
}
