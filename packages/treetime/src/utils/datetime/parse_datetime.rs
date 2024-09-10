use crate::make_error;
use crate::utils::datetime::datetime::{date_from_iso, date_from_rfc2822};
use crate::utils::datetime::options::DateParserOptions;
use chrono::{DateTime, NaiveDateTime, Utc};
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
    .map(|naive_datetime| DateTime::<Utc>::from_utc(naive_datetime, Utc))
    .wrap_err_with(|| format!("When parsing datetime '{datetime_str}' using format '{format}'"))
}

const DATETIME_FORMATS: &[&str] = &[
  // ISO 8601 datetime (2024-07-23T18:37:20)
  // Well-known as the international standard format used widely in web protocols.
  "%Y-%m-%dT%H:%M:%S",
  // ISO 8601 with fractional seconds (2024-07-23 18:37:20.215)
  // Used in contexts requiring high precision time stamping.
  "%Y-%m-%d %H:%M:%S%.f",
  // ISO 8601 with timezone offset (2024-07-23 18:37:20 +0000)
  // Essential for international applications dealing with multiple time zones.
  "%Y-%m-%d %H:%M:%S %z",
  // ISO 8601 with colon in timezone offset (2024-07-23T18:37:20.215Z)
  // Standard format in many web and network protocols.
  "%Y-%m-%dT%H:%M:%S%:z",
  // POSIX locale time format (Tue Jul 23 18:37:20 2024)
  // Used in POSIX systems for locale-aware applications.
  "%a %b %e %H:%M:%S %Y",
  // ctime() format (Tue Jul 23 18:37:20 2024)
  // Common format in UNIX and similar systems.
  "%a %b %d %H:%M:%S %Y",
  // European datetime with seconds (23-07-2024 18:37:20)
  // Common format in European countries.
  "%d-%m-%Y %H:%M:%S",
  // Day abbreviated month year time (23 Jul 2024 18:37:20)
  // Commonly seen in informal communications and logs.
  "%d %b %Y %H:%M:%S",
  // Compact numerical datetime (20240723183720)
  // Often used in compact data formats and logging.
  "%Y%m%d%H%M%S",
  // Year month day hour minute second, dot-separated (202407231837.20)
  // Uncommon, may be used in specialized data processing.
  "%Y%m%d%H%M.%S",
  // Year month day hour minute second with timezone offset, no separators (20240723183720 +0000)
  // Useful in systems that require compact timezone-aware timestamps.
  "%Y%m%d%H%M%S %z",
  // Year month day hour minute second with UTC offset, no separators (20240723183720UTC+0000)
  // Specific use in systems requiring explicit UTC designation.
  "%Y%m%d%H%MUTC%z",
  // Year month day hour minute second with Zulu offset, no separators (20240723183720Z+0000)
  // Common in military and aviation contexts for Zulu (GMT) time notation.
  "%Y%m%d%H%MZ%z",
  // Year month day hour minute second with timezone offset, hyphenated (2024-07-23 18:37 +0000)
  // Variant format for readability with explicit timezone.
  "%Y-%m-%d %H:%M %z",
  // Dot-separated European date and time (23.07.2024 18:37:20)
  // Common in many Eastern and Central European countries.
  "%d.%m.%Y %H:%M:%S",
  // US datetime with 12-hour clock (07-23-2024 06:37 PM)
  // Typically used in the United States.
  "%m-%d-%Y %I:%M %p",
];
