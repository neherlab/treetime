use crate::{make_error, make_report};
use chrono::{DateTime, Datelike, FixedOffset, NaiveDate, NaiveTime, Utc};
use eyre::{Report, WrapErr};

pub(crate) const NOON: NaiveTime = NaiveTime::from_hms_opt(12, 0, 0).expect("12:00:00 is a valid time");
pub(crate) const LAST_NANOSECOND: NaiveTime =
  NaiveTime::from_hms_nano_opt(23, 59, 59, 999_999_999).expect("23:59:59.999999999 is a valid time");
const MILLISECONDS_PER_SECOND: i64 = 1000;

#[expect(
  clippy::unwrap_used,
  reason = "parses date literals written in code, where a malformed literal is a programming error"
)]
pub fn iso(date_str: impl AsRef<str>) -> DateTime<Utc> {
  date_from_iso(date_str).unwrap()
}

pub fn date_from_rfc2822(date_str: impl AsRef<str>) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  let parsed = DateTime::<FixedOffset>::parse_from_rfc2822(date_str)
    .wrap_err_with(|| format!("When parsing datetime '{date_str}' using RFC 2822 format"))?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

pub fn date_from_iso(date_str: impl AsRef<str>) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  let parsed = DateTime::<FixedOffset>::parse_from_rfc3339(date_str)
    .wrap_err_with(|| format!("When parsing datetime '{date_str}' using RFC 3339 (ISO 8601) format"))?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

pub fn date_now() -> DateTime<Utc> {
  Utc::now()
}

pub fn date_to_timestamp(datetime: &DateTime<Utc>) -> i64 {
  datetime.timestamp() * MILLISECONDS_PER_SECOND
}

pub fn date_format(datetime: &DateTime<Utc>) -> String {
  datetime.format("%Y-%m-%d %H:%M:%S").to_string()
}

pub fn date_format_precise(datetime: &DateTime<Utc>) -> String {
  datetime.format("%Y-%m-%d %H:%M:%S%.3f").to_string()
}

pub fn timestamp_to_date(timestamp: i64) -> Result<DateTime<Utc>, Report> {
  DateTime::from_timestamp_millis(timestamp)
    .ok_or_else(|| make_report!("Timestamp {timestamp} ms is outside the supported date range"))
}

pub fn date_format_safe(datetime: &DateTime<Utc>) -> String {
  datetime.format("%Y-%m-%d_%H-%M-%S").to_string()
}

pub fn days_in_month(year: u32, month: u32) -> Result<u32, Report> {
  if !(1..=12).contains(&month) {
    return make_error!("Invalid month: {month}");
  }
  let (next_year, next_month) = if month == 12 { (year + 1, 1) } else { (year, month + 1) };
  let signed_next_year =
    i32::try_from(next_year).wrap_err_with(|| format!("Year {next_year} is outside the supported date range"))?;
  let last_day = NaiveDate::from_ymd_opt(signed_next_year, next_month, 1)
    .and_then(|first_of_next_month| first_of_next_month.pred_opt())
    .ok_or_else(|| make_report!("Year {year} is outside the supported date range"))?;
  Ok(last_day.day())
}
