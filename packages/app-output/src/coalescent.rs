use eyre::Report;
use std::path::Path;
use treetime::timetree::coalescent::CoalescentOutput;
use treetime_io::csv::CsvStructFileWriter;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn write_coalescent_json(output: &CoalescentOutput, path: impl AsRef<Path>) -> Result<(), Report> {
  json_write_file(path, output, JsonPretty(true))
}

pub fn write_coalescent_delimited(
  output: &CoalescentOutput,
  path: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut writer = CsvStructFileWriter::new(path, delimiter)?;
  for row in output.rows() {
    writer.write(&row)?;
  }
  writer.finish()
}
