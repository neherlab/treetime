use crate::commands::homoplasy::args::TreetimeHomoplasyArgs;
use crate::make_error;
use eyre::Report;

pub fn run_homoplasy(_: TreetimeHomoplasyArgs) -> Result<(), Report> {
  make_error!("The homoplasy command is not yet implemented")
}
