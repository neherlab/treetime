use crate::csv::{detect_csv_delimiter, get_col_name, normalize_csv_headers};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::{make_internal_report, make_report};

pub fn read_discrete_attrs_from_reader<T>(
  reader: impl Read,
  delimiter: u8,
  name_candidates: &[String],
  name_column: &Option<String>,
  value_column: &Option<String>,
  parser: impl Fn(&str) -> Result<T, Report>,
) -> Result<(BTreeMap<String, T>, String), Report> {
  let mut reader = CsvReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .comment(None)
    .from_reader(reader);

  let headers = reader
    .headers()
    .map(normalize_csv_headers)
    .map_err(|err| make_report!("{err}"))?;

  let name_column_idx = get_col_name(&headers, name_candidates, name_column)?;
  let value_column_idx = get_col_name(&headers, &[], value_column)?;

  let value_name = headers[value_column_idx].clone();

  let values = reader
    .records()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      convert_record::<T>(index, &record, name_column_idx, value_column_idx, &parser)
    })
    .collect::<Result<BTreeMap<String, T>, Report>>()?;

  Ok((values, value_name))
}

pub fn read_discrete_attrs_from_str<T>(
  content: &str,
  delimiter: u8,
  name_candidates: &[String],
  name_column: &Option<String>,
  value_column: &Option<String>,
  parser: impl Fn(&str) -> Result<T, Report>,
) -> Result<(BTreeMap<String, T>, String), Report> {
  read_discrete_attrs_from_reader(
    content.as_bytes(),
    delimiter,
    name_candidates,
    name_column,
    value_column,
    parser,
  )
}

pub fn read_discrete_attrs<T>(
  filepath: impl AsRef<Path>,
  delimiters: &[char],
  name_candidates: &[String],
  name_column: &Option<String>,
  value_column: &Option<String>,
  parser: impl Fn(&str) -> Result<T, Report>,
) -> Result<(BTreeMap<String, T>, String), Report> {
  let filepath = filepath.as_ref();
  let mut file =
    open_file_or_stdin(&Some(filepath)).wrap_err_with(|| format!("When reading file: '{}'", filepath.display()))?;
  let delimiter = detect_csv_delimiter(&mut *file, filepath, delimiters, |headers| {
    get_col_name(headers, name_candidates, name_column).is_ok() && get_col_name(headers, &[], value_column).is_ok()
  })
  .wrap_err_with(|| format!("When detecting CSV delimiter for '{}'", filepath.display()))?;
  read_discrete_attrs_from_reader(file, delimiter, name_candidates, name_column, value_column, parser)
    .wrap_err_with(|| format!("When reading discrete attributes from file: '{}'", filepath.display()))
}

pub fn convert_record<T>(
  index: usize,
  record: &StringRecord,
  name_column_idx: usize,
  value_column_idx: usize,
  parser: &impl Fn(&str) -> Result<T, Report>,
) -> Result<(String, T), Report> {
  let key = record
    .get(name_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{name_column_idx}'"))?
    .to_owned();

  let value = record
    .get(value_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{value_column_idx}'"))?;

  let value = parser(value)?;

  Ok((key, value))
}
