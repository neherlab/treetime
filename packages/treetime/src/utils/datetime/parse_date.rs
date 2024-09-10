use crate::make_error;
use crate::utils::datetime::date_range::DateRange;
use crate::utils::datetime::options::{DateParserOptions, TimeOfDay};
use crate::utils::datetime::year_frac::year_fraction_to_date;
use chrono::{DateTime, NaiveDate, Utc};
use eyre::{Report, WrapErr};
use lazy_static::lazy_static;
use regex::Regex;

pub fn parse_date(date_str: impl AsRef<str>, options: &DateParserOptions) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();

  if let Ok(date) = parse_date_with_formats(date_str, DATE_FORMATS, options) {
    return Ok(date);
  }

  if let Ok(year_fraction) = date_str.parse::<f64>() {
    if year_fraction.is_finite() {
      return Ok(year_fraction_to_date(year_fraction));
    }
  }

  make_error!("Unrecognized date format: {date_str}")
}

pub fn parse_date_with_formats(
  date_str: impl AsRef<str>,
  date_formats: impl IntoIterator<Item = impl AsRef<str>>,
  options: &DateParserOptions,
) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  for format in date_formats {
    let format = format.as_ref();
    if let Ok(datetime) = parse_date_with_format(date_str, format, options) {
      return Ok(datetime);
    }
  }
  make_error!("Unrecognized date format: {date_str}")
}

pub fn parse_date_with_format(
  date_str: impl AsRef<str>,
  format: impl AsRef<str>,
  options: &DateParserOptions,
) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  let format = format.as_ref();
  NaiveDate::parse_from_str(date_str, format)
    .map(|naive_date| match options.default_time_of_day {
      TimeOfDay::Dawn => naive_date.and_hms(0, 0, 0),
      TimeOfDay::Noon => naive_date.and_hms(12, 0, 0),
      TimeOfDay::Dusk => naive_date.and_hms_nano(23, 59, 59, 999_999_999),
      TimeOfDay::Custom(time) => naive_date.and_time(time),
      TimeOfDay::CustomFn(func) => naive_date.and_time(func(&naive_date)),
    })
    .map(|naive_datetime| DateTime::<Utc>::from_utc(naive_datetime, Utc))
    .wrap_err_with(|| format!("When parsing date '{date_str}' using format '{format}'"))
}

pub const DATE_FORMATS: &[&str] = &[
  // Y, m, d
  "%Y-%m-%d",
  "%Y/%m/%d",
  "%Y.%m.%d",
  // d, m, Y
  "%d-%m-%Y",
  "%d/%m/%Y",
  "%d.%m.%Y",
  //
  "%Y%m%d",
  // Year and day of the year (2024-205)
  "%Y-%j",
  // ISO week date (2024-W30-2)
  "%Y-W%W-%w",
  // ISO week and week year (30-2024)
  "%V-%G",
  // US date format (07/23/2024)
  "%m/%d/%Y",
  // Month name day, year (July 23, 2024)
  "%B %d, %Y",
];

pub fn parse_date_range(date_range_str: &str, options: &DateParserOptions) -> Result<DateRange, Report> {
  for (regex, fmt) in DATE_RANGE_REGEX.iter() {
    if let Some(captures) = regex.captures(date_range_str) {
      if let (Some(begin), Some(end)) = (captures.name("begin"), captures.name("end")) {
        let begin = parse_date_with_format(begin.as_str(), fmt, options);
        let end = parse_date_with_format(end.as_str(), fmt, options);
        if let (Ok(begin), Ok(end)) = (begin, end) {
          return Ok(DateRange::new(begin, end));
        }
      }
    }
  }
  make_error!("Unrecognized date range format: {date_range_str}")
}

lazy_static! {
  static ref DATE_RANGE_REGEX: Vec<(Regex, String)> = create_date_range_regexes();
}

fn create_date_range_regexes() -> Vec<(Regex, String)> {
  #[rustfmt::skip]
  const DATE_PATTERNS: &[(&str, &str)] = &[
    ( r"\d{4}-\d{2}-\d{2}",     "%Y-%m-%d" ),
    ( r"\d{4}\d{2}\d{2}",       "%Y%m%d"   ),
    ( r"\d{2}\.\d{2}\.\d{4}",   "%d.%m.%Y" ),
    ( r"\d{4}-W\d{2}-7",        "%Y-W%W-7" ),
  ];

  let mut regexes = vec![];

  // Separator-delimited formats, e.g. "2021-12-17 / 2024-07-23"
  {
    const SEPARATORS: &[&str] = &["/", "..", "...", "-", ",", ";", ":"];
    for pattern in DATE_PATTERNS {
      for separator in SEPARATORS {
        let format = pattern.0;
        let regex = Regex::new(&format!(r"^(?P<begin>{format})\s*{separator}\s*(?P<end>{format})$")).unwrap();
        regexes.push((regex, pattern.1.to_owned()));
      }
    }
  }

  // Parenthesized formats, e.g. "[ 2021-12-17 , 2024-07-23 ]"
  {
    const PARENS: &[(&str, &str)] = &[(r"\[", r"\]"), (r"\(", r"\)"), (r"\{", r"\}")];
    const SEPARATORS: &[&str] = &["..", "...", ",", ";", ":"];
    for pattern in DATE_PATTERNS {
      for (left_paren, right_paren) in PARENS {
        for separator in SEPARATORS {
          let format = pattern.0;
          let regex = Regex::new(&format!(
            r"^{left_paren}\s*(P<begin>{format})\s*{right_paren}\s*{separator}\s*{left_paren}\s*(P<end>{format})\s*{right_paren}$"
          ))
          .unwrap();
          regexes.push((regex, pattern.1.to_owned()));
        }
      }
    }
  }

  regexes
}
