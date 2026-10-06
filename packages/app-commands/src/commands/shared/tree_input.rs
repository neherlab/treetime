use crate::commands::shared::input_warnings::warn_duplicate_names;
use eyre::Report;
use std::path::Path;
use treetime::progress::{LogSink, RunWarningKind};
use treetime_io::nwk::NwkParse;
use treetime_io::tree::tree_read_file;

pub fn read_input_tree(path: &Path, log: &dyn LogSink) -> Result<NwkParse, Report> {
  let parse = tree_read_file(path)?;
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateNodeNames,
    &format!(
      "The tree '{}' gives the same name to more than one node:",
      path.display()
    ),
    "Nodes with the same name receive the same data from the other inputs, and augur node data keeps one entry per name.",
    &parse.duplicate_names,
  );
  Ok(parse)
}
