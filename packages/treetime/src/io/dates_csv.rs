use crate::io::csv::{get_col_name, guess_csv_delimiter};
use crate::io::file::open_file_or_stdin;
use crate::utils::datetime::{date_from_formats, date_from_iso, date_from_rfc2822, date_to_year_fraction};
use crate::{make_error, make_internal_report, vec_of_owned};
use csv::{ReaderBuilder as CsvReaderBuilder, StringRecord, Trim};
use eyre::{eyre, Report, WrapErr};
use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::LazyLock;

#[derive(Clone, Copy, Debug)]
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
  let file = open_file_or_stdin(&Some(filepath)).wrap_err_with(|| format!("When reading file: {filepath:#?}"))?;
  let delimiter =
    guess_csv_delimiter(filepath).wrap_err_with(|| format!("When guessing CSV delimiter for {filepath:#?}"))?;

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
    .wrap_err_with(|| format!("When detecting name column in {filepath:#?}"))?;

  let date_column_idx = get_col_name(&headers, &vec_of_owned!["date"], date_column)
    .wrap_err_with(|| format!("When detecting date column in {filepath:#?}"))?;

  reader
    .records()
    .enumerate()
    .map(|(index, record)| {
      let record = record?;
      convert_record(index, &record, name_column_idx, date_column_idx)
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

  let date = parse_date(date).wrap_err_with(|| format!("Row '{index}': When parsing date column"))?;

  Ok((name, date))
}

pub fn parse_date(date_str: &str) -> Result<Option<DateOrRange>, Report> {
  let date_str = &date_str.trim_matches(|c: char| !c.is_ascii() || c.is_whitespace());

  if is_nil(date_str) {
    Ok(None)
  }
  // Try to parse as year fraction: `2022.7`
  else if let Ok(date) = date_str.parse::<f64>() {
    if date.is_finite() {
      Ok(Some(DateOrRange::YearFraction(date)))
    } else {
      Ok(None)
    }
  }
  // Try to parse as year fraction range: `[2022.6:2022.7]`
  else if let Ok(range) = parse_date_range(date_str) {
    Ok(Some(DateOrRange::YearFractionRange(range)))
  }
  // Try to parse as ISO 8601 (RFC 3339) date: `2022-07-22T16:40:59Z`
  else if let Ok(date) = date_from_iso(date_str) {
    Ok(Some(DateOrRange::YearFraction(date_to_year_fraction(&date))))
  }
  // Try to parse as RFC 2822 date: `Wed, 22 Jul 2022 18:48:09 GMT`
  else if let Ok(date) = date_from_rfc2822(date_str) {
    Ok(Some(DateOrRange::YearFraction(date_to_year_fraction(&date))))
  }
  // Try to parse as various date formats
  else if let Ok(date) = date_from_formats(date_str) {
    Ok(Some(DateOrRange::YearFraction(date_to_year_fraction(&date))))
  } else {
    // Give up
    make_error!("Unknown date format: {date_str}")
  }
}

fn is_nil(input: &str) -> bool {
  const NON_VALUES: &[&str] = &[
    "na",
    "n/a",
    "nan",
    "null",
    "nil",
    "none",
    "empty",
    "missing",
    "undefined",
  ];
  static REGEX: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\?+|-+)$").unwrap());
  let input_lower = input.to_lowercase();
  input.is_empty() || NON_VALUES.contains(&input_lower.as_str()) || REGEX.is_match(&input_lower)
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

#[cfg(test)]
mod tests {
  #![allow(clippy::excessive_precision)]
  use super::*;
  use crate::pretty_assert_ulps_eq;
  use chrono::{DateTime, TimeZone, Utc};
  use rstest::*;

  fn ymd(year: i32, month: u32, day: u32) -> DateTime<Utc> {
    Utc.ymd(year, month, day).and_hms(0, 0, 0)
  }

  #[rustfmt::skip]
  #[rstest]
  // TODO: implement BC dates? (double-check test expectations)
  // #[case::start_of_mayan_long_count_calendar          (ymd(-3114,  8, 11), -3113.0027397260)]
  // #[case::eclipse_of_ho_and_hi                        (ymd(-2137, 10, 22), -2136.8095890410)]
  // #[case::beginning_of_reign_of_thutmose_iii          (ymd(-1479,  4, 28), -1424.1931506849)]
  // #[case::founding_of_rome                            (ymd( -753,  4, 21),  -752.2246575342)]
  // #[case::birth_of_confucius                          (ymd( -551,  9, 28),  -550.7965753424)]
  // #[case::assassination_of_julius_caesar              (ymd(  -44,  3, 15),   -43.2054794520)]
  // #[case::end_of_bc                                   (ymd(   -1, 12, 31),     0.0013698630)]
  #[case::start_of_ad                                 (ymd(    1,  1,  1),     1.0013698630)]
  #[case::pompeii_destruction                         (ymd(   79,  8, 24),    79.6452054794)]
  #[case::founding_of_baghdad                         (ymd(  762,  7, 30),   762.5767123287)]
  #[case::supernova_1054_ad                           (ymd( 1054,  7,  4),  1054.5054794520)]
  #[case::first_crusade                               (ymd( 1096,  8, 15),  1096.6215846994)]
  #[case::gregorian_calendar_introduced               (ymd( 1582, 10, 15),  1582.7876712328)]
  #[case::completion_of_the_suez_canal                (ymd( 1869, 11, 17),  1869.8780821917)]
  #[case::wright_brothers_first_flight                (ymd( 1903, 12, 17),  1903.9602739726)]
  #[case::discovery_of_penicillin                     (ymd( 1928, 9,  28),  1928.7418032786)]
  #[case::invention_of_the_transistor                 (ymd( 1947, 12, 23),  1947.9767123287)]
  #[case::intel_microprocessor_release                (ymd( 1971,  9, 6),   1971.6808219178)]
  #[case::first_iss_launch                            (ymd( 1998, 11, 20),  1998.8863013698)]
  #[case::dot_com_bubble_burst                        (ymd( 2000,  3, 10),  2000.1898907103)]
  #[case::iphone_launch                               (ymd( 2007,  6, 29),  2007.4917808219)]
  #[case::this_test_written                           (ymd( 2024,  9,  5),  2024.6789617486)]
  #[case::halley_comet_return                         (ymd( 2061,  7, 28),  2061.5712328767)]
  #[case::transit_of_venus                            (ymd( 2117, 12, 11),  2117.9438356164)]
  #[case::pluto_full_orbit_since_discovery            (ymd( 2178,  3, 23),  2178.2232876712)]
  #[case::thousand_years_since_this_test              (ymd( 3024,  9,  5),  3024.6789617486)]
  //
  #[case::half_century                                (ymd( 1950,  7,  2),  1950.5000000000)]
  #[case::quarter_century                             (ymd( 2025,  4,  1),  2025.2479452054)]
  #[case::three_quarters_century                      (ymd( 2075,  9, 30),  2075.7465753424)]
  #[case::hapl_millennium                             (ymd( 2500, 12, 31),  2500.9986301369)]
  #[case::halfway_non_leap_year                       (ymd( 2023,  7,  2),  2023.5000000000)]
  //
  #[case::leap_year_2024_day_before                   (ymd( 2024,  2, 28),  2024.15983606553)]
  #[case::leap_year_2024_leap_day                     (ymd( 2024,  2, 29),  2024.16256830600)]
  #[case::leap_year_2024_day_after                    (ymd( 2024,  3,  1),  2024.16530054648)]
  #[case::leap_year_2024_halfway                      (ymd( 2024,  7,  2),  2024.50136612028)]
  #[case::leap_year_2024_last_day                     (ymd( 2024, 12, 31),  2024.99863387971)]
  #[case::leap_year_2024_day_after_year               (ymd( 2025,  1,  1),  2025.00136986303)]
  //
  #[case::round_leap_year_2000_day_before             (ymd( 2000,  2, 28),  2000.15983606557)]
  #[case::round_leap_year_2000_leap_day               (ymd( 2000,  2, 29),  2000.16256830601)]
  #[case::round_leap_year_2000_day_after              (ymd( 2000,  3,  1),  2000.16530054644)]
  #[case::round_leap_year_2000_halfway                (ymd( 2000,  7,  2),  2000.50136612021)]
  #[case::round_leap_year_2000_last_day               (ymd( 2000, 12, 31),  2000.99863387978)]
  #[case::round_leap_year_2000_day_after_year         (ymd( 2001,  1,  1),  2001.00136986301)]
  //
  #[case::excluded_leap_year_1900_day_before          (ymd( 1900,  2, 27),  1900.15753424657)]
  #[case::excluded_leap_year_1900_want_to_be_leap_day (ymd( 1900,  2, 28),  1900.16027397260)]
  #[case::excluded_leap_year_1900_day_after           (ymd( 1900,  3,  1),  1900.16301369863)]
  #[case::excluded_leap_year_1900_last_day            (ymd( 1900, 12, 31),  1900.99863013698)]
  #[case::excluded_leap_year_1900_halfway             (ymd( 1900,  7,  2),  1900.50000000000)]
  #[case::excluded_leap_year_1900_day_after_year      (ymd( 1901,  1,  1),  1901.00136986301)]
  //
  #[case::non_leap_year_2023_day_before               (ymd( 2023,  2, 27),  2023.15753424657)]
  #[case::non_leap_year_2023_want_to_be_leap_day      (ymd( 2023,  2, 28),  2023.16027397260)]
  #[case::non_leap_year_2023_day_after                (ymd( 2023,  3,  1),  2023.16301369863)]
  #[case::non_leap_year_2023_halfway                  (ymd( 2023,  7,  2),  2023.50000000000)]
  #[case::non_leap_year_2023_last_day                 (ymd( 2023, 12, 31),  2023.99863013698)]
  #[case::non_leap_year_2023_day_after_year           (ymd( 2024,  1,  1),  2024.00136612021)]
  //
  #[case::y2k_day_before                              (ymd( 1999, 12, 31),  1999.99863013698)]
  #[case::y2k                                         (ymd( 2000,  1,  1),  2000.00136612021)]
  #[case::y2k_day_after                               (ymd( 2000,  1,  2),  2000.00409836065)]
  #[case::y2k38_before                                (ymd( 2038,  1, 17),  2038.04520547945)]
  #[case::y2k38                                       (ymd( 2038,  1, 18),  2038.04794520547)]
  #[case::y2k38_after                                 (ymd( 2038,  1, 19),  2038.05068493150)]
  //
  #[case::unix_epoch_day_before                       (ymd( 1969, 12, 31),  1969.99863013698)]
  #[case::unix_epoch                                  (ymd( 1970,  1,  1),  1970.00136986301)]
  #[case::unix_epoch_day_after                        (ymd( 1970,  1,  2),  1970.00410958904)]
  //
  #[trace]
  fn test_date_to_year_fraction(#[case] input: DateTime<Utc>, #[case] expected: f64) {
    let actual = date_to_year_fraction(&input);
    pretty_assert_ulps_eq!(expected, actual, epsilon = 1e-6);
  }
}

// #[rstest]
// #[case("2000-03-29", "%Y-%m-%d", 2000.242)]
// fn test_get_numerical_date_from_value_not_ambiguous(
//   #[case] input: &str,
//   #[case] format: &str,
//   #[case] expected: f64,
// ) {
//   assert_eq!(dates::get_numerical_date_from_value(input, format), expected);
// }
//
// #[rstest]
// #[case("2000-01-XX", "%Y-%m-%d", 2000.001, 2000.083)]
// fn test_get_numerical_date_from_value_ambiguous_day(
//   #[case] input: &str,
//   #[case] format: &str,
//   #[case] min_expected: f64,
//   #[case] max_expected: f64,
// ) {
//   let (min_date, max_date) = dates::get_numerical_date_from_value(input, format);
//   assert_eq!(min_date, min_expected);
//   assert_eq!(max_date, max_expected);
// }
//
// #[rstest]
// #[case("2000-XX-XX", "%Y-%m-%d", 2000.001, 2000.999)]
// fn test_get_numerical_date_from_value_ambiguous_month_and_day(
//   #[case] input: &str,
//   #[case] format: &str,
//   #[case] min_expected: f64,
//   #[case] max_expected: f64,
// ) {
//   let (min_date, max_date) = dates::get_numerical_date_from_value(input, format);
//   assert_eq!(min_date, min_expected);
//   assert_eq!(max_date, max_expected);
// }
//
// #[rstest]
// #[case("200X", "min", 2000)]
// #[case("200X", "max", 2009)]
// #[case("20X0", "max", 2090)]
// #[case("X000", "max", 9000)]
// #[case("XXXX", "min", 1)]
// #[case("XXXX", "max", 9999)]
// fn test_resolve_uncertain_int(#[case] input: &str, #[case] min_or_max: &str, #[case] expected: u32) {
//   let result = resolve_uncertain_int(input, min_or_max).unwrap();
//   assert_eq!(result, expected);
// }
//
// #[rstest]
// #[case("2000-01-01", (NaiveDate::from_ymd(2000, 1, 1), NaiveDate::from_ymd(2000, 1, 1)))]
// #[case("2000-02-XX", (NaiveDate::from_ymd(2000, 2, 1), NaiveDate::from_ymd(2000, 2, 29)))]
// #[case("2000-XX-XX", (NaiveDate::from_ymd(2000, 1, 1), NaiveDate::from_ymd(2000, 12, 31)))]
// fn test_range(#[case] input: &str, #[case] expected_range: (NaiveDate, NaiveDate)) {
//   let ad = AmbiguousDate::new(input, "%Y-%m-%d");
//   assert_eq!(ad.range(None).unwrap(), expected_range);
// }
//
// #[rstest]
// #[case("2000-01-01", "%Y-%m-%d", NaiveDate::from_ymd(2000, 1, 1))]
// #[case("2000-01-XX", "%Y-%m-%d", NaiveDate::from_ymd(2000, 1, 1))]
// #[case("2000-XX-XX", "%Y-%m-%d", NaiveDate::from_ymd(2000, 1, 1))]
// fn test_uncertain_date_components(#[case] input: &str, #[case] format: &str, #[case] expected_date: NaiveDate) {
//   let ad = AmbiguousDate::new(input, format);
//   let components = ad.uncertain_date_components().unwrap();
//   assert_eq!(
//     NaiveDate::from_ymd(
//       components["Y"].parse::<i32>().unwrap(),
//       components["m"].parse::<u32>().unwrap(),
//       components["d"].parse::<u32>().unwrap()
//     ),
//     expected_date
//   );
// }
//
// #[test]
// fn test_uncertain_date_components_error() {
//   let result = AmbiguousDate::new("5-5-5-5-5", "%Y-%m-%d").uncertain_date_components();
//   assert!(result.is_err());
// }
//
// #[rstest]
// #[case("2000-03-XX", "%Y-%m-%d", (NaiveDate::from_ymd(2000, 3, 1), NaiveDate::from_ymd(2000, 3, 31)))]
// #[case("2020-XX-XX", None, (NaiveDate::from_ymd(2020, 1, 1), NaiveDate::from_ymd(2020, 12, 31)))]
// #[case("2020-12-XX", None, (NaiveDate::from_ymd(2020, 12, 1), NaiveDate::from_ymd(2020, 12, 31)))]
// #[case("2020-12-01", None, (NaiveDate::from_ymd(2020, 12, 1), NaiveDate::from_ymd(2020, 12, 1)))]
// #[case("20XX-XX-XX", Some((2010, None)), (NaiveDate::from_ymd(2010, 1, 1), NaiveDate::from_ymd(2010, 12, 31)))]
// fn test_range_with_bounds(
//   #[case] input: &str,
//   #[case] min_max_year: Option<(i32, Option<i32>)>,
//   #[case] expected_range: (NaiveDate, NaiveDate),
// ) {
//   let ad = AmbiguousDate::new(input, "%Y-%m-%d");
//   assert_eq!(ad.range(min_max_year).unwrap(), expected_range);
// }
//
// #[rstest]
// #[case("200X-01-01", "Year contains uncertainty, so month and day must also be uncertain")]
// #[case("2000-XX-01", "Month contains uncertainty, so day must also be uncertain")]
// fn test_assert_only_less_significant_uncertainty(#[case] input: &str, #[case] expected_error: &str) {
//   let ad = AmbiguousDate::new(input, "%Y-%m-%d");
//   let result = ad.uncertain_date_components();
//   assert!(result.is_err());
//   assert_eq!(format!("{}", result.err().unwrap()), expected_error);
// }
//
// #[rstest]
// #[case("20XX-XX-XX", [1950, 1960], "Not possible for date to fall within bounds")]
// #[case("19XX-XX-XX", [2000, 2010], "Not possible for date to fall within bounds")]
// fn test_min_max_year_date_error(#[case] input: &str, #[case] min_max_year: [i32; 2], #[case] expected_error: &str) {
//   let ad = AmbiguousDate::new(input, "%Y-%m-%d");
//   let result = ad.range(Some((min_max_year[0], Some(min_max_year[1]))));
//   assert!(result.is_err());
//   assert_eq!(format!("{}", result.err().unwrap()), expected_error);
// }

// use chrono::{Datelike, NaiveDate};
// use std::collections::HashMap;
// use std::fmt;
//
// #[derive(Debug)]
// pub enum DateError {
//   #[error("invalid date: `{0}`")]
//   InvalidDate(String),
//
//   #[error("year bounds invalid: `{0}`")]
//   InvalidYearBounds(String),
// }
//
// fn tuple_to_date(year: i32, month: u32, day: u32) -> NaiveDate {
//   let month = month.min(12);
//   let day = day.min(max_day_for_year_month(year, month));
//
//   NaiveDate::from_ymd(year, month, day)
// }
//
// fn max_day_for_year_month(year: i32, month: u32) -> u32 {
//   NaiveDate::from_ymd(year, month, 1)
//     .with_month0(month - 1)
//     .unwrap()
//     .last_day()
// }
//
// fn resolve_uncertain_int(uncertain_string: &str, min_or_max: &str) -> Result<u32, DateError> {
//   let result = match min_or_max {
//     "min" => uncertain_string.replace('X', "0"),
//     "max" => uncertain_string.replace('X', "9"),
//     _ => {
//       return Err(DateError::InvalidDate(format!(
//         "Invalid min_or_max parameter: {}",
//         min_or_max
//       )))
//     }
//   };
//
//   let int_result = result.parse::<u32>().or_else(|_| {
//     Err(DateError::InvalidDate(format!(
//       "Cannot parse integer from: {}",
//       uncertain_string
//     )))
//   })?;
//
//   Ok(int_result.max(1))
// }
//
// struct AmbiguousDate {
//   uncertain_date: String,
//   fmt: String,
// }
//
// impl AmbiguousDate {
//   pub fn new(uncertain_date: &str, fmt: &str) -> Self {
//     Self {
//       uncertain_date: uncertain_date.to_string(),
//       fmt: fmt.to_string(),
//     }
//   }
//
//   pub fn range(&self, min_max_year: Option<(i32, Option<i32>)>) -> Result<(NaiveDate, NaiveDate), DateError> {
//     let components = self.uncertain_date_components()?;
//     let min_date = tuple_to_date(
//       resolve_uncertain_int(components.get("Y").unwrap(), "min")? as i32,
//       resolve_uncertain_int(components.get("m").unwrap(), "min")?,
//       resolve_uncertain_int(components.get("d").unwrap(), "min")?,
//     );
//
//     let max_date = tuple_to_date(
//       resolve_uncertain_int(components.get("Y").unwrap(), "max")? as i32,
//       resolve_uncertain_int(components.get("m").unwrap(), "max")?,
//       resolve_uncertain_int(components.get("d").unwrap(), "max")?,
//     );
//
//     let today = chrono::Local::today().naive_local();
//
//     let (min_date, max_date) = match min_max_year {
//       Some((lower, upper)) => {
//         let lower_bound = NaiveDate::from_ymd(lower, 1, 1);
//         let upper_bound = upper.map_or_else(|| today, |year| NaiveDate::from_ymd(year, 12, 31));
//
//         (min_date.max(lower_bound), max_date.min(upper_bound))
//       }
//       None => (min_date, max_date),
//     };
//
//     Ok((min_date.min(today), max_date.min(today)))
//   }
//
//   fn uncertain_date_components(&self) -> Result<HashMap<&str, &str>, DateError> {
//     let re = Regex::new(&self.regex()?)
//       .map_err(|_| DateError::InvalidDate(format!("Invalid regex for format: {}", self.fmt)))?;
//     let captures = re
//       .captures(&self.uncertain_date)
//       .ok_or_else(|| DateError::InvalidDate(format!("Date does not match format `{}`.", self.fmt)))?;
//
//     let mut components = HashMap::new();
//     for component in self.fmt_components().iter() {
//       components.insert(*component, captures.name(component).unwrap().as_str());
//     }
//
//     Ok(components)
//   }
//
//   fn fmt_components(&self) -> Vec<&str> {
//     self.fmt.split('%').filter(|&s| !s.is_empty()).collect()
//   }
//
//   fn regex(&self) -> Result<String, DateError> {
//     Ok(
//       self
//         .fmt
//         .replace("%Y", "(?P<Y>....)")
//         .replace("%m", "(?P<m>..?)")
//         .replace("%d", "(?P<d>..?)"),
//     )
//   }
//}
