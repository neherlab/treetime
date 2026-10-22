use crate::output_format::table_format;
use deser::Serialize;
use deser::de::DeserializeOwned;
use eyre::Report;
use std::path::Path;
use treetime_io::csv::{CsvWriter, TableFormat, csv_read_file, csv_write_file};
use treetime_utils::io::file::FileWriter;
use treetime_utils::make_internal_report;

pub fn table_write_file<T: Serialize>(filepath: &Path, rows: impl IntoIterator<Item = T>) -> Result<(), Report> {
  csv_write_file(filepath, rows, planned_table_format(filepath)?)
}

pub fn table_create(filepath: &Path) -> Result<CsvWriter<FileWriter>, Report> {
  CsvWriter::create(filepath, planned_table_format(filepath)?)
}

pub fn table_read_file<T: DeserializeOwned>(filepath: &Path) -> Result<Vec<T>, Report> {
  csv_read_file(filepath, planned_table_format(filepath)?)
}

fn planned_table_format(filepath: &Path) -> Result<TableFormat, Report> {
  table_format(filepath)
    .ok_or_else(|| make_internal_report!("No table format for '{}' after output planning", filepath.display()))
}
