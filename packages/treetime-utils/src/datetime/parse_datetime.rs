use crate::datetime::datetime::{date_from_iso, date_from_rfc2822};
use crate::datetime::options::DateParserOptions;
use crate::make_error;
use chrono::{DateTime, NaiveDateTime, TimeZone, Utc};
use eyre::{Context, Report};

pub fn parse_datetime(datetime_str: impl AsRef<str>, options: &DateParserOptions) -> Result<DateTime<Utc>, Report> {
  parse_datetime_with_formats(datetime_str, DATETIME_FORMATS, options)
}

pub fn parse_datetime_with_formats(
  datetime_str: impl AsRef<str>,
  datetime_formats: impl IntoIterator<Item = impl AsRef<str>>,
  _options: &DateParserOptions,
) -> Result<DateTime<Utc>, Report> {
  let date_str = datetime_str.as_ref();

  if let Ok(dt) = date_from_iso(date_str) {
    return Ok(dt.with_timezone(&Utc));
  }

  if let Ok(dt) = date_from_rfc2822(date_str) {
    return Ok(dt.with_timezone(&Utc));
  }

  for format in datetime_formats {
    let format = format.as_ref();
    if let Ok(datetime) = parse_datetime_with_format(date_str, format) {
      return Ok(datetime);
    }
  }
  make_error!("Unrecognized date format: {date_str}")
}

pub fn parse_datetime_with_format(
  datetime_str: impl AsRef<str>,
  format: impl AsRef<str>,
) -> Result<DateTime<Utc>, Report> {
  let datetime_str = datetime_str.as_ref();
  let format = format.as_ref();
  NaiveDateTime::parse_from_str(datetime_str, format)
    .map(|naive_datetime| Utc.from_utc_datetime(&naive_datetime))
    .wrap_err_with(|| format!("When parsing datetime '{datetime_str}' using format '{format}'"))
}

const DATETIME_FORMATS: &[&str] = &[
  "%Y-%m-%dT%H:%M:%S",
  "%Y-%m-%d %H:%M:%S%.f",
  "%Y-%m-%d %H:%M:%S %z",
  "%Y-%m-%dT%H:%M:%S%:z",
  "%a %b %e %H:%M:%S %Y",
  "%a %b %d %H:%M:%S %Y",
  "%d-%m-%Y %H:%M:%S",
  "%d %b %Y %H:%M:%S",
  "%Y%m%d%H%M%S",
  "%Y%m%d%H%M.%S",
  "%Y%m%d%H%M%S %z",
  "%Y%m%d%H%MUTC%z",
  "%Y%m%d%H%MZ%z",
  "%Y-%m-%d %H:%M %z",
  "%d.%m.%Y %H:%M:%S",
  "%m-%d-%Y %I:%M %p",
];
