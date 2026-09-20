use crate::make_error;
use chrono::{DateTime, Datelike, FixedOffset, NaiveDate, TimeZone, Utc};
use eyre::{Report, WrapErr};

pub fn date_now() -> DateTime<Utc> {
  Utc::now()
}

// Parse RFC 3339 and ISO 8601 datetime string (1996-12-19T16:39:57-08:00)
pub fn date_from_iso(date_str: impl AsRef<str>) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  let parsed = DateTime::<FixedOffset>::parse_from_rfc3339(date_str)
    .wrap_err_with(|| format!("When parsing datetime '{date_str}' using RFC 3339 (ISO 8601) format"))?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

#[allow(clippy::unwrap_used, reason = "unwrap on a value an upstream invariant guarantees is present")]
pub fn iso(date_str: impl AsRef<str>) -> DateTime<Utc> {
  date_from_iso(date_str).unwrap()
}

// Parse RFC 2822 datetime string (Tue, 1 Jul 2003 10:52:37 +0200)
pub fn date_from_rfc2822(date_str: impl AsRef<str>) -> Result<DateTime<Utc>, Report> {
  let date_str = date_str.as_ref();
  let parsed = DateTime::<FixedOffset>::parse_from_rfc2822(date_str)
    .wrap_err_with(|| format!("When parsing datetime '{date_str}' using RFC 2822 format"))?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

pub fn date_to_timestamp(datetime: &DateTime<Utc>) -> i64 {
  datetime.timestamp() * 1000
}

#[allow(clippy::expect_used, reason = "expect on a value an upstream invariant guarantees is present")]
/// Convert millisecond timestamp to DateTime.
///
/// # Panics
/// Panics if timestamp is out of range for DateTime (approximately ±262,000 years from epoch).
pub fn timestamp_to_date(timestamp: i64) -> DateTime<Utc> {
  DateTime::from_timestamp_millis(timestamp).expect("timestamp out of DateTime range")
}

pub fn timestamp_from_iso(iso: &str) -> Result<i64, Report> {
  let datetime = date_from_iso(iso)?;
  let timestamp = date_to_timestamp(&datetime);
  Ok(timestamp)
}

pub fn timestamp_now() -> i64 {
  date_to_timestamp(&date_now())
}

pub fn date_format(datetime: &DateTime<Utc>) -> String {
  datetime.format("%Y-%m-%d %H:%M:%S").to_string()
}

pub fn timestamp_format(timestamp: i64) -> String {
  date_format(&timestamp_to_date(timestamp))
}

pub fn date_format_precise(datetime: &DateTime<Utc>) -> String {
  datetime.format("%Y-%m-%d %H:%M:%S%.3f").to_string()
}

pub fn timestamp_format_precise(timestamp: i64) -> String {
  date_format_precise(&timestamp_to_date(timestamp))
}

pub fn date_format_safe(datetime: &DateTime<Utc>) -> String {
  datetime.format("%Y-%m-%d_%H-%M-%S").to_string()
}

pub fn timestamp_format_safe(timestamp: i64) -> String {
  date_format_safe(&timestamp_to_date(timestamp))
}

pub fn ymd(year: i32, month: u32, day: u32) -> DateTime<Utc> {
  Utc.with_ymd_and_hms(year, month, day, 0, 0, 0).unwrap()
}

#[allow(clippy::as_conversions, clippy::expect_used, reason = "count/index numeric cast is exact for the domain range; expect on a value an upstream invariant guarantees is present")]
pub fn days_in_month(year: u32, month: u32) -> Result<u32, Report> {
  if !(1..=12).contains(&month) {
    return make_error!("Invalid month: {month}");
  }
  let (next_year, next_month) = if month == 12 { (year + 1, 1) } else { (year, month + 1) };
  let last_day = NaiveDate::from_ymd_opt(next_year as i32, next_month, 1)
    .and_then(|first_of_next_month| first_of_next_month.pred_opt())
    .expect("a validated 1-12 month always has a last day");
  Ok(last_day.day())
}
