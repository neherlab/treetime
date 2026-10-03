use crate::auspice_types::AuspiceTree;
use eyre::{Report, WrapErr};
use std::io::Write;
use std::path::Path;
use treetime_utils::io::file::write_file_with;
use treetime_utils::io::json::{JsonPretty, json_write};

pub fn auspice_write_file(filepath: impl AsRef<Path>, tree: &AuspiceTree) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  write_file_with(filepath, |f| {
    auspice_write(&mut *f, tree)?;
    writeln!(f).map_err(Report::new)
  })
  .wrap_err_with(|| format!("When writing Auspice v2 JSON file '{}'", filepath.display()))
}

fn auspice_write(writer: &mut impl Write, tree: &AuspiceTree) -> Result<(), Report> {
  json_write(writer, tree, JsonPretty(true)).wrap_err("When writing Auspice v2 JSON")
}
