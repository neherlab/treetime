use crate::auspice_types::AuspiceTree;
use eyre::{Report, WrapErr};
use std::io::Write;
use std::path::Path;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::json::{JsonPretty, json_write};

pub fn auspice_write_file(filepath: impl AsRef<Path>, tree: &AuspiceTree) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let context = || format!("When writing Auspice v2 JSON file '{}'", filepath.display());
  let mut f = create_file_or_stdout(filepath)?;
  auspice_write(&mut f, tree).wrap_err_with(context)?;
  writeln!(f).wrap_err_with(context)
}

fn auspice_write(writer: &mut impl Write, tree: &AuspiceTree) -> Result<(), Report> {
  json_write(writer, tree, JsonPretty(true)).wrap_err("When writing Auspice v2 JSON")
}
