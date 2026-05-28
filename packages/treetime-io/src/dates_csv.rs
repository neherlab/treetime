use crate::csv::{get_col_name, guess_csv_delimiter};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr};
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
use treetime_utils::{make_internal_report, make_report, vec_of_owned};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DateExact {
  pub value: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DateRange {
  pub start: f64,
  pub end: f64,
}

impl DateRange {
  pub fn width(&self) -> f64 {
    self.end - self.start
  }

  pub fn contains(&self, value: f64) -> bool {
    value >= self.start && value <= self.end
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum DateValue {
  Exact(DateExact),
  Uncertain(DateRange),
  Range(DateRange),
}

impl DateValue {
  #[inline]
  pub fn mean(&self) -> f64 {
    match self {
      DateValue::Exact(d) => d.value,
      DateValue::Uncertain(r) | DateValue::Range(r) => f64::midpoint(r.start, r.end),
    }
  }

  pub fn is_exact(&self) -> bool {
    matches!(self, DateValue::Exact(_))
  }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DateConstraint {
  pub raw: String,
  pub value: DateValue,
}

impl DateConstraint {
  pub fn exact(value: f64) -> Self {
    Self {
      raw: value.to_string(),
      value: DateValue::Exact(DateExact { value }),
    }
  }

  #[inline]
  pub fn mean(&self) -> f64 {
    self.value.mean()
  }

  pub fn is_exact(&self) -> bool {
    self.value.is_exact()
  }
}

pub type DatesMap = BTreeMap<String, Option<DateConstraint>>;
pub type DateRecord = (String, Option<DateConstraint>);

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
    .map_err(|err| make_report!("{err}"))?
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

pub fn read_date(date_str: &str, options: &DateParserOptions) -> Result<Option<DateConstraint>, Report> {
  let trimmed = date_str.trim();

  if let Ok(year_fraction) = trimmed.parse::<f64>()
    && year_fraction.is_finite()
    && (1000.0..=3000.0).contains(&year_fraction)
  {
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Exact(DateExact { value: year_fraction }),
    }));
  }

  if let Ok(date) = parse_date(trimmed, options) {
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Exact(DateExact {
        value: date_to_year_fraction(&date),
      }),
    }));
  }

  if let Ok(date_range) = parse_date_uncertain(trimmed, options) {
    let (start, end) = date_range_to_year_fraction_range(&date_range);
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Uncertain(DateRange { start, end }),
    }));
  }

  if let Ok(date_range) = parse_date_range(trimmed, options) {
    let (start, end) = date_range_to_year_fraction_range(&date_range);
    return Ok(Some(DateConstraint {
      raw: trimmed.to_owned(),
      value: DateValue::Range(DateRange { start, end }),
    }));
  }

  Ok(None)
}
