use crate::csv::{
  ColumnRule, delimiter_from_path, detect_csv_delimiter, find_column, normalize_csv_headers, table_records,
};
use eyre::{Report, WrapErr};
use std::io::BufRead;
use std::path::Path;
use treetime_utils::io::file::read_file_with;
use treetime_utils::make_internal_report;

pub fn discrete_attrs_read_file<T>(
  filepath: impl AsRef<Path>,
  delimiters: &[char],
  name_candidates: &[String],
  name_column: Option<&str>,
  value_column: Option<&str>,
  parser: impl Fn(&str) -> Result<T, Report>,
) -> Result<(Vec<(String, T)>, String), Report> {
  let filepath = filepath.as_ref();
  read_file_with(filepath, |reader| {
    discrete_attrs_read(
      reader,
      delimiter_from_path(filepath),
      delimiters,
      name_candidates,
      name_column,
      value_column,
      parser,
    )
  })
}

pub fn discrete_attrs_read<T>(
  mut reader: impl BufRead,
  path_delimiter: Option<u8>,
  delimiters: &[char],
  name_candidates: &[String],
  name_column: Option<&str>,
  value_column: Option<&str>,
  parser: impl Fn(&str) -> Result<T, Report>,
) -> Result<(Vec<(String, T)>, String), Report> {
  let name_rule = ColumnRule::requested_or(name_column, ColumnRule::OneOf(name_candidates));
  let value_rule = ColumnRule::requested_or(value_column, ColumnRule::OneOf(&[]));
  let delimiter = detect_csv_delimiter(&mut reader, path_delimiter, delimiters, |headers| {
    find_column(headers, name_rule).is_ok() && find_column(headers, value_rule).is_ok()
  })
  .wrap_err("When detecting the table delimiter")?;
  let mut reader = table_records(reader, delimiter);

  let written = reader.read::<Vec<String>>()?.unwrap_or_default();
  let headers = normalize_csv_headers(&written);

  let name_column_idx = find_column(&written, name_rule)?;
  let value_column_idx = find_column(&written, value_rule)?;

  let value_name = headers[value_column_idx].clone();

  let values = reader
    .iter::<Vec<String>>()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      convert_record::<T>(index, &record, name_column_idx, value_column_idx, &parser)
    })
    .collect::<Result<Vec<(String, T)>, Report>>()?;

  Ok((values, value_name))
}

fn convert_record<T>(
  index: usize,
  record: &[String],
  name_column_idx: usize,
  value_column_idx: usize,
  parser: &impl Fn(&str) -> Result<T, Report>,
) -> Result<(String, T), Report> {
  let key = record
    .get(name_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{name_column_idx}'"))?
    .clone();

  let value = record
    .get(value_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{value_column_idx}'"))?;

  let value = parser(value)?;

  Ok((key, value))
}
