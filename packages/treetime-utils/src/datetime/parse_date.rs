use crate::datetime::date_range::DateRange;
use crate::datetime::options::{DateParserOptions, TimeOfDay};
use crate::datetime::year_fraction::year_fraction_to_date;
use crate::make_error;
use chrono::{DateTime, NaiveDate, TimeZone, Utc};
use eyre::{Report, WrapErr};
use regex::Regex;
use std::sync::LazyLock;

pub fn parse_date(date_str: impl AsRef<str>, options: &DateParserOptions) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();

  if let Ok(date) = parse_date_with_formats(date_str, DATE_FORMATS, options) {
    return Ok(date);
  }

  if let Ok(year_fraction) = date_str.parse::<f64>()
    && year_fraction.is_finite()
  {
    return Ok(year_fraction_to_date(year_fraction));
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

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
pub fn parse_date_with_format(
  date_str: impl AsRef<str>,
  format: impl AsRef<str>,
  options: &DateParserOptions,
) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  let format = format.as_ref();
  NaiveDate::parse_from_str(date_str, format)
    .map(|naive_date| match options.default_time_of_day {
      TimeOfDay::Dawn => naive_date.and_hms_opt(0, 0, 0).unwrap(),
      TimeOfDay::Noon => naive_date.and_hms_opt(12, 0, 0).unwrap(),
      TimeOfDay::Dusk => naive_date.and_hms_nano_opt(23, 59, 59, 999_999_999).unwrap(),
      TimeOfDay::Custom(time) => naive_date.and_time(time),
      TimeOfDay::CustomFn(func) => naive_date.and_time(func(&naive_date)),
    })
    .map(|naive_datetime| Utc.from_utc_datetime(&naive_datetime))
    .wrap_err_with(|| format!("When parsing date '{date_str}' using format '{format}'"))
}

pub const DATE_FORMATS: &[&str] = &[
  "%Y-%m-%d",
  "%Y/%m/%d",
  "%Y.%m.%d",
  "%d-%m-%Y",
  "%d/%m/%Y",
  "%d.%m.%Y",
  "%Y%m%d",
  "%Y-%j",
  "%Y-W%W-%w",
  "%V-%G",
  "%m/%d/%Y",
  "%B %d, %Y",
];

pub fn parse_date_range(date_range_str: &str, options: &DateParserOptions) -> Result<DateRange, Report> {
  for (regex, fmt) in DATE_RANGE_REGEX.iter() {
    if let Some(captures) = regex.captures(date_range_str)
      && let (Some(begin), Some(end)) = (captures.name("begin"), captures.name("end"))
    {
      let begin = parse_date_with_format(begin.as_str(), fmt, options);
      let end = parse_date_with_format(end.as_str(), fmt, options);
      if let (Ok(begin), Ok(end)) = (begin, end) {
        return Ok(DateRange::new(begin, end));
      }
    }
  }
  make_error!("Unrecognized date range format: {date_range_str}")
}

static DATE_RANGE_REGEX: LazyLock<Vec<(Regex, String)>> = LazyLock::new(create_date_range_regexes);

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn create_date_range_regexes() -> Vec<(Regex, String)> {
  #[rustfmt::skip]
  const DATE_PATTERNS: &[(&str, &str)] = &[
    ( r"\d{4}-\d{2}-\d{2}",     "%Y-%m-%d" ),
    ( r"\d{4}\d{2}\d{2}",       "%Y%m%d"   ),
    ( r"\d{2}\.\d{2}\.\d{4}",   "%d.%m.%Y" ),
    ( r"\d{4}-W\d{2}-7",        "%Y-W%W-7" ),
  ];

  let mut regexes = vec![];

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

  {
    const PARENS: &[(&str, &str)] = &[(r"\[", r"\]"), (r"\(", r"\)"), (r"\{", r"\}")];
    const SEPARATORS: &[&str] = &["..", "...", ",", ";", ":"];
    for pattern in DATE_PATTERNS {
      for (left_paren, right_paren) in PARENS {
        for separator in SEPARATORS {
          let format = pattern.0;
          let regex = Regex::new(&format!(
            r"^{left_paren}\s*(?P<begin>{format})\s*{right_paren}\s*{separator}\s*{left_paren}\s*(?P<end>{format})\s*{right_paren}$"
          ))
          .unwrap();
          regexes.push((regex, pattern.1.to_owned()));
        }
      }
    }
  }

  regexes
}
