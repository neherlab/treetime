use eyre::Report;
use std::io::Write;
use std::path::Path;
use treetime::timetree::coalescent::CoalescentOutput;
use treetime_io::csv::CsvStructWriter;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn write_coalescent_json(output: &CoalescentOutput, path: impl AsRef<Path>) -> Result<(), Report> {
  json_write_file(path, output, JsonPretty(true))
}

pub fn write_coalescent_delimited(
  output: &CoalescentOutput,
  path: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let file = create_file_or_stdout(path)?;
  write_coalescent_delimited_to(output, file, delimiter)?.finish()
}

pub fn write_coalescent_delimited_to<W: Write + Send>(
  output: &CoalescentOutput,
  writer: W,
  delimiter: u8,
) -> Result<W, Report> {
  let mut writer = CsvStructWriter::new(writer, delimiter)?;
  for row in output.rows() {
    writer.write(&row)?;
  }
  writer.into_inner()
}
