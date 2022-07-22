use crate::io::fs::extension;
use crate::utils::datetime::{date_from_formats, date_from_iso, date_from_rfc2822, date_to_year_fraction};
use crate::{make_error, make_internal_report, vec_of_owned};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{eyre, Report, WrapErr};
use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Clone, Debug)]
pub enum DateOrRange {
  YearFraction(f64),
  YearFractionRange((f64, f64)),
}

#[derive(Clone, Debug)]
pub struct DateRecord {
  name: String,
  date: DateOrRange,
}

pub fn read_dates(
  filepath: impl AsRef<Path>,
  name_column: &Option<String>,
  date_column: &Option<String>,
) -> Result<Vec<DateRecord>, Report> {
  let filepath = filepath.as_ref();
  let file = BufReader::new(File::open(&filepath)?);
  let delimiter =
    guess_csv_delimiter(&filepath).wrap_err_with(|| format!("When guessing CSV delimiter for '{filepath:#?}'"))?;

  let mut reader = CsvReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .from_reader(file);

  let headers = reader
    .headers()
    .map(|record| record.iter().map(str::to_owned).collect_vec())
    .map_err(|err| eyre!("{err}"))?;

  let name_column_idx = get_col_name(&headers, &vec_of_owned!["name", "strain", "accession"], name_column)
    .wrap_err_with(|| format!("When detecting name column in '{filepath:#?}'"))?;

  let date_column_idx = get_col_name(&headers, &vec_of_owned!["date"], date_column)
    .wrap_err_with(|| format!("When detecting date column in '{filepath:#?}'"))?;

  reader
    .records()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      convert_record(index, &record, name_column_idx, date_column_idx)
    })
    .collect::<Result<Vec<DateRecord>, Report>>()
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
    .ok_or_else(|| make_internal_report!("Row '{index}': Unable to get column with index '{name_column_idx}'"))?;

  let date = parse_date(date).wrap_err_with(|| format!("Row '{index}': When parsing date column"))?;

  Ok(DateRecord { name, date })
}

pub fn parse_date(date_str: &str) -> Result<DateOrRange, Report> {
  // Try to parse as year fraction: `2022.7`
  if let Ok(date) = date_str.parse::<f64>() {
    Ok(DateOrRange::YearFraction(date))
  }
  // Try to parse as year fraction range: `[2022.6:2022.7]`
  else if let Ok(range) = parse_date_range(date_str) {
    Ok(DateOrRange::YearFractionRange(range))
  }
  // Try to parse as ISO 8601 (RFC 3339) date: `2022-07-22T16:40:59Z`
  else if let Ok(date) = date_from_iso(date_str) {
    Ok(DateOrRange::YearFraction(date_to_year_fraction(&date)))
  }
  // Try to parse as RFC 2822 date: `Wed, 22 Jul 2022 18:48:09 GMT`
  else if let Ok(date) = date_from_rfc2822(date_str) {
    Ok(DateOrRange::YearFraction(date_to_year_fraction(&date)))
  }
  // Try to parse as various date formats
  else if let Ok(date) = date_from_formats(date_str) {
    Ok(DateOrRange::YearFraction(date_to_year_fraction(&date)))
  } else {
    // Give up
    make_error!("Unknown date format: {date_str}")
  }
}

pub fn parse_date_range(date_str: &str) -> Result<(f64, f64), Report> {
  const REGEX_STR: &str = r"\[(?P<begin>\d{4}\.\d+):(?P<end>\d{4}\.\d+)]";
  lazy_static! {
    static ref RE: Regex = Regex::new(REGEX_STR)
      .wrap_err_with(|| format!("When compiling regular expression '{REGEX_STR}'"))
      .unwrap();
  }

  if let Some(captures) = RE.captures(date_str) {
    match (captures.name("begin"), captures.name("end")) {
      (Some(begin), Some(end)) => {
        let begin = begin.as_str().parse::<f64>()?;
        let end = end.as_str().parse::<f64>()?;
        Ok((begin, end))
      }
      _ => make_error!("Unable to parse date range: '{date_str}'"),
    }
  } else {
    make_error!("Unable to parse date range: '{date_str}'")
  }
}

pub fn get_col_name(
  headers: &[String],
  possible_names: &[String],
  provided_name: &Option<String>,
) -> Result<usize, Report> {
  if let Some(provided_name) = provided_name {
    match headers.iter().position(|header| header == provided_name) {
      Some(idx) => Ok(idx),
      None => make_error!("Requested column '{provided_name}' not found"),
    }
  } else {
    headers
      .iter()
      .position(|header| possible_names.contains(header))
      .ok_or_else(|| {
        eyre!(
          "Unable to find column:\n  Looking for: {}\n  Found: {}",
          possible_names.join(", "),
          headers.join(",")
        )
      })
  }
}

pub fn guess_csv_delimiter(filepath: impl AsRef<Path>) -> Result<u8, Report> {
  let filepath = filepath.as_ref();
  let ext = extension(filepath)
    .ok_or_else(|| eyre!("Unable to detect file extension: '{filepath:?}': "))?
    .to_lowercase();
  match ext.as_str() {
    "csv" => Ok(b','),
    "tsv" => Ok(b'\t'),
    "ssv" => Ok(b';'),
    _ => make_error!("Unknown file extension: '{ext}'"),
  }
}
