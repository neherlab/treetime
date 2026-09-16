use treetime::timetree::coalescent::CoalescentOutput;
use eyre::{Report, WrapErr};
use std::path::Path;
use treetime_io::csv::{CsvStructFileWriter, CsvStructWriter};
use treetime_utils::io::json::{JsonPretty, json_write_file, json_write_str};
use treetime_utils::make_report;

/// Serialize the coalescent output as one rich JSON document (`inputs` + `outputs.segments`).
pub fn coalescent_json_str(output: &CoalescentOutput) -> Result<String, Report> {
  json_write_str(output, JsonPretty(true))
}

/// Serialize the coalescent segments as flat delimited rows (header + one row per segment).
///
/// The header is derived from the serde field names of [`CoalescentSegmentRow`], so the dotted
/// names (`segment.start`, `T_c.value`, ...) become the column titles. `delimiter` selects CSV
/// (`b','`) or TSV (`b'\t'`).
pub fn coalescent_delimited_str(output: &CoalescentOutput, delimiter: u8) -> Result<String, Report> {
  let mut writer = CsvStructWriter::new(Vec::<u8>::new(), delimiter)?;
  for row in output.rows() {
    writer.write(&row)?;
  }
  let bytes = writer
    .writer
    .into_inner()
    .map_err(|err| make_report!("Coalescent delimited serialization failed to flush: {err}"))?;
  String::from_utf8(bytes).wrap_err("Coalescent delimited output is not valid UTF-8")
}

/// Write the rich coalescent document to a JSON file.
pub fn write_coalescent_json(output: &CoalescentOutput, path: impl AsRef<Path>) -> Result<(), Report> {
  json_write_file(path, output, JsonPretty(true))
}

/// Write the flat coalescent segments to a delimited file (CSV with `b','`, TSV with `b'\t'`).
pub fn write_coalescent_delimited(
  output: &CoalescentOutput,
  path: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut writer = CsvStructFileWriter::new(path, delimiter)?;
  for row in output.rows() {
    writer.write(&row)?;
  }
  Ok(())
}
