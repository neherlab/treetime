use chrono::{DateTime, FixedOffset, TimeZone, Utc};
use eyre::{Report, WrapErr};
use std::time::{Duration, UNIX_EPOCH};
use time::util::days_in_year_month;
use time::Month;

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

pub fn timestamp_to_date(timestamp: i64) -> DateTime<Utc> {
  DateTime::<Utc>::from(UNIX_EPOCH + Duration::from_millis(timestamp as u64))
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
  Utc.ymd(year, month, day).and_hms(0, 0, 0)
}

pub fn days_in_month(year: u32, month: u32) -> Result<u32, Report> {
  let month_obj = Month::try_from(month as u8)?;
  let days = days_in_year_month(year as i32, month_obj) as u32;
  Ok(days)
}
