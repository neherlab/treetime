use crate::io::csv::{get_col_name, guess_csv_delimiter};
use crate::io::file::open_file_or_stdin;
use crate::utils::datetime::options::DateParserOptions;
use crate::utils::datetime::parse_date::{parse_date, parse_date_range};
use crate::utils::datetime::parse_uncertain_date::parse_date_uncertain;
use crate::utils::datetime::year_frac::{date_range_to_year_fraction_range, date_to_year_fraction};
use crate::{make_internal_report, vec_of_owned};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{Report, WrapErr, eyre};
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq)]
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

  let mut reader = CsvReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
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

  let date_column_idx = get_col_name(&headers, &vec_of_owned!["date"], date_column)
    .wrap_err_with(|| format!("When detecting date column in '{}'", filepath.display()))?;

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
  if let Ok(date) = parse_date(date_str, options) {
    Ok(Some(DateOrRange::YearFraction(date_to_year_fraction(&date))))
  } else if let Ok(date_range) = parse_date_uncertain(date_str, options) {
    let yf_range = date_range_to_year_fraction_range(&date_range);
    Ok(Some(DateOrRange::YearFractionRange(yf_range)))
  } else if let Ok(date_range) = parse_date_range(date_str, options) {
    let yf_range = date_range_to_year_fraction_range(&date_range);
    Ok(Some(DateOrRange::YearFractionRange(yf_range)))
  } else {
    Ok(None)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case("")]
  #[case("   ")]
  #[case("NaN")]
  #[case("null")]
  #[case("XXXX-XX-XX")]
  #[case("XXXX/XX/XX")]
  #[case("XXXX.XX.XX")]
  #[trace]
  fn test_date_read_empty(#[case] input: &str) -> Result<(), Report>  {
    let actual = read_date(input, &DateParserOptions::default())?;
    assert_eq!(None, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("2024-07-23",   2024.5587431693)]
  #[case("2024/07/23",   2024.5587431693)]
  #[case("2024.07.23",   2024.5587431693)]
  //
  #[case("2024-07-XX",   2024.5396174863)]
  #[case("2024/07/XX",   2024.5396174863)]
  #[case("2024.07.XX",   2024.5396174863)]
  //
  #[case("2024-XX-XX",   2024.5000000000)]
  #[case("2024/XX/XX",   2024.5000000000)]
  #[case("2024.XX.XX",   2024.5000000000)]
  //
  #[case("2024.558743",  2024.5587431693)]
  #[case("20240723",     2024.5587431693)]
  #[trace]
  fn test_date_read_ok(#[case] input: &str, #[case] expected: f64) -> Result<(), Report>  {
    let actual = read_date(input, &DateParserOptions::default())?.unwrap().mean();
    pretty_assert_ulps_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }
}
