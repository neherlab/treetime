use crate::commands::homoplasy::args::TreetimeHomoplasyArgs;
use crate::commands::homoplasy::result::HomoplasyResult;
use crate::make_error;
use eyre::Report;

pub fn run_homoplasy(_: TreetimeHomoplasyArgs) -> Result<HomoplasyResult, Report> {
  make_error!("The homoplasy command is not yet implemented in v1")
}
