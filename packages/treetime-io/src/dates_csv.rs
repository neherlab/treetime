use crate::csv::{detect_csv_delimiter, get_col_name, normalize_csv_headers};
use csv::{ReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr};
use std::io::Read;
use std::path::Path;
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::{parse_date, parse_date_range};
use treetime_utils::datetime::parse_uncertain_date::parse_date_uncertain;
use treetime_utils::datetime::year_fraction::{date_range_to_year_fraction_range, date_to_year_fraction};
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::{make_internal_report, make_report, vec_of_owned};

pub use treetime_primitives::date::{DateConstraint, DateExact, DateRange, DateValue, DatesMap};

pub(crate) fn read_dates_from_str(
  content: &str,
  delimiter: u8,
  name_candidates: &[String],
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let reader = content.as_bytes();
  read_dates_from_reader(
    reader,
    delimiter,
    name_candidates,
    name_column.as_deref(),
    date_column.as_deref(),
  )
}

pub fn read_dates(
  filepath: impl AsRef<Path>,
  delimiters: &[char],
  name_candidates: &[String],
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<DatesMap, Report> {
  let filepath = filepath.as_ref();
  let mut file =
    open_file_or_stdin(&Some(filepath)).wrap_err_with(|| format!("When reading file: '{}'", filepath.display()))?;
  let delimiter = detect_csv_delimiter(&mut *file, filepath, delimiters, |headers| {
    get_col_name(headers, name_candidates, name_column.as_deref()).is_ok()
      && get_col_name(headers, &vec_of_owned!["date"], date_column.as_deref()).is_ok()
  })
  .wrap_err_with(|| format!("When detecting CSV delimiter for '{}'", filepath.display()))?;
  read_dates_from_reader(
    file,
    delimiter,
    name_candidates,
    name_column.as_deref(),
    date_column.as_deref(),
  )
  .wrap_err_with(|| format!("When reading dates from file: '{}'", filepath.display()))
}

fn read_dates_from_reader(
  reader: impl Read,
  delimiter: u8,
  name_candidates: &[String],
  name_column: Option<&str>,
  date_column: Option<&str>,
) -> Result<DatesMap, Report> {
  let mut reader = ReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .from_reader(reader);

  let headers = reader
    .headers()
    .map(normalize_csv_headers)
    .map_err(|err| make_report!("{err}"))?;

  let name_column_idx = get_col_name(&headers, name_candidates, name_column)?;
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

fn convert_record(
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

type DateRecord = (String, Option<DateConstraint>);

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    error_dropped_by_pattern,
    reason = "each date notation is tried in turn; a value that matches none is reported by the caller as missing"
  )
)]
pub(crate) fn read_date(date_str: &str, options: &DateParserOptions) -> Result<Option<DateConstraint>, Report> {
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
