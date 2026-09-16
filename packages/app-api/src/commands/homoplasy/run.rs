use crate::commands::homoplasy::args::TreetimeHomoplasyArgs;
use crate::commands::homoplasy::result::HomoplasyResult;
use eyre::Report;
use treetime::make_error;

pub fn run_homoplasy(_: TreetimeHomoplasyArgs) -> Result<HomoplasyResult, Report> {
  make_error!("The homoplasy command is not yet implemented in v1")
}
