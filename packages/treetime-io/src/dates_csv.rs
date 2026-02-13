use crate::csv::{get_col_name, guess_csv_delimiter};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr, eyre};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::{parse_date, parse_date_range};
use treetime_utils::datetime::parse_uncertain_date::parse_date_uncertain;
use treetime_utils::datetime::year_fraction::{date_range_to_year_fraction_range, date_to_year_fraction};
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::{make_internal_report, vec_of_owned};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum DateOrRange {
  YearFraction(f64),
  YearFractionRange((f64, f64)),
}

impl DateOrRange {
  #[inline]
  pub fn mean(&self) -> f64 {
    match self {
      DateOrRange::YearFraction(date) => *date,
      DateOrRange::YearFractionRange((begin, end)) => (begin + end) / 2.0,
    }
  }
}

pub type DatesMap = BTreeMap<String, Option<DateOrRange>>;
pub type DateRecord = (String, Option<DateOrRange>);

pub fn read_dates_from_reader(
  reader: impl Read,
  delimiter: u8,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let mut reader = CsvReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .from_reader(reader);

  let headers = reader
    .headers()
    .map(|record| record.iter().map(str::to_owned).collect_vec())
    .map_err(|err| eyre!("{err}"))?
    .iter()
    .map(|header| header.trim_start_matches('#').trim_end_matches('#').trim().to_owned())
    .collect_vec();

  let name_column_idx = get_col_name(&headers, &vec_of_owned!["name", "strain", "accession"], name_column)?;
  let date_column_idx = get_col_name(&headers, &vec_of_owned!["date"], date_column)?;

  reader
    .records()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      convert_record(index, &record, name_column_idx, date_column_idx)
        .wrap_err_with(|| format!("When reading row {index}, column '{date_column_idx}'"))
    })
    .collect::<Result<DatesMap, Report>>()
}

pub fn read_dates_from_str(
  content: &str,
  delimiter: u8,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let reader = content.as_bytes();
  read_dates_from_reader(reader, delimiter, name_column, date_column)
}

pub fn read_dates(
  filepath: impl AsRef<Path>,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let filepath = filepath.as_ref();
  let file =
    open_file_or_stdin(&Some(filepath)).wrap_err_with(|| format!("When reading file: '{}'", filepath.display()))?;
  let delimiter = guess_csv_delimiter(filepath)
    .wrap_err_with(|| format!("When guessing CSV delimiter for '{}'", filepath.display()))?;
  read_dates_from_reader(file, delimiter, name_column, date_column)
    .wrap_err_with(|| format!("When reading dates from file: '{}'", filepath.display()))
}

pub fn convert_record(
  index: usize,
  record: &StringRecord,
  name_column_idx: usize,
  date_column_idx: usize,
) -> Result<DateRecord, Report> {
  let name = record
    .get(name_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{name_column_idx}'"))?
    .to_owned();

  let date = record
    .get(date_column_idx)
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{date_column_idx}'"))?;

  let date = read_date(date, &DateParserOptions::default())?;

  Ok((name, date))
}

pub fn read_date(date_str: &str, options: &DateParserOptions) -> Result<Option<DateOrRange>, Report> {
  let trimmed = date_str.trim();

  // Try parsing as year fraction first (e.g. "2020.5", "2003.84052019")
  // Only accept values in reasonable year range (1000-3000) to avoid parsing compact dates like "20240723" as floats
  if let Ok(year_fraction) = trimmed.parse::<f64>()
    && year_fraction.is_finite()
    && (1000.0..=3000.0).contains(&year_fraction)
  {
    return Ok(Some(DateOrRange::YearFraction(year_fraction)));
  }

  // Try parsing as formatted date (e.g. "2020-07-15", "20200715", "2020/07/15")
  if let Ok(date) = parse_date(trimmed, options) {
    return Ok(Some(DateOrRange::YearFraction(date_to_year_fraction(&date))));
  }

  // Try parsing as uncertain date (e.g. "2020-XX-XX")
  if let Ok(date_range) = parse_date_uncertain(trimmed, options) {
    return Ok(Some(DateOrRange::YearFractionRange(date_range_to_year_fraction_range(
      &date_range,
    ))));
  }

  // Try parsing as date range (e.g. "2020-01-01/2020-12-31")
  if let Ok(date_range) = parse_date_range(trimmed, options) {
    return Ok(Some(DateOrRange::YearFractionRange(date_range_to_year_fraction_range(
      &date_range,
    ))));
  }

  // Unable to parse
  Ok(None)
}
