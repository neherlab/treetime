use crate::io::csv::{get_col_name, guess_csv_delimiter};
use crate::io::file::open_file_or_stdin;
use crate::{make_internal_report, vec_of_owned};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr, eyre};
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::Path;

pub fn read_discrete_attrs<T>(
  filepath: impl AsRef<Path>,
  name_column: &Option<String>,
  value_column: &Option<String>,
  parser: impl Fn(&str) -> Result<T, Report>,
) -> Result<(BTreeMap<String, T>, String), Report> {
  let filepath = filepath.as_ref();
  let file =
    open_file_or_stdin(&Some(filepath)).wrap_err_with(|| format!("When reading file: '{}'", filepath.display()))?;
  let delimiter = guess_csv_delimiter(filepath)
    .wrap_err_with(|| format!("When guessing CSV delimiter for '{}'", filepath.display()))?;

  let mut reader = CsvReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .comment(None)
    .from_reader(file);

  let headers = reader
    .headers()
    .map(|record| record.iter().map(str::to_owned).collect_vec())
    .map_err(|err| eyre!("{err}"))?
    .iter()
    .map(|header| header.trim_start_matches('#').trim_end_matches('#').to_owned())
    .collect_vec();

  let name_column_idx = get_col_name(&headers, &vec_of_owned!["name", "strain", "accession"], name_column)
    .wrap_err_with(|| format!("When detecting name column in '{}'", filepath.display()))?;

  let value_column_idx = get_col_name(&headers, &[], value_column)
    .wrap_err_with(|| format!("When detecting attribute column in '{}'", filepath.display()))?;

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
