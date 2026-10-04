use crate::output_plan::OutputSelection;
use eyre::Report;
use serde::Serialize;
use serde::de::DeserializeOwned;
use std::path::Path;
use treetime_io::csv::{CsvWriter, TableFormat, csv_read_file, csv_write_file};
use treetime_utils::io::file::FileWriter;
use treetime_utils::make_internal_report;

pub fn table_write_file<T: Serialize>(
  selection: OutputSelection,
  filepath: &Path,
  rows: impl IntoIterator<Item = T>,
) -> Result<(), Report> {
  csv_write_file(filepath, rows, table_format(selection)?)
}

pub fn table_create(selection: OutputSelection, filepath: &Path) -> Result<CsvWriter<FileWriter>, Report> {
  CsvWriter::create(filepath, table_format(selection)?)
}

pub fn table_read_file<T: DeserializeOwned>(selection: OutputSelection, filepath: &Path) -> Result<Vec<T>, Report> {
  csv_read_file(filepath, table_format(selection)?)
}

fn table_format(selection: OutputSelection) -> Result<TableFormat, Report> {
  selection
    .table_format()
    .ok_or_else(|| make_internal_report!("Output '{}' is not a table", selection.as_ref()))
}
