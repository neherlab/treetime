use crate::make_error;
use chrono::{DateTime, Datelike, FixedOffset, NaiveDate, Utc};
use eyre::Report;
use itertools::Itertools;
use std::time::{Duration, UNIX_EPOCH};
use time::util::days_in_year;

pub fn date_now() -> DateTime<Utc> {
  Utc::now()
}

pub fn date_from_iso(iso: &str) -> Result<DateTime<Utc>, Report> {
  let parsed = DateTime::<FixedOffset>::parse_from_rfc3339(iso)?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

pub fn date_from_rfc2822(date_str: &str) -> Result<DateTime<Utc>, Report> {
  let parsed = DateTime::<FixedOffset>::parse_from_rfc2822(date_str)?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

pub fn date_from_format(format: &str, date_str: &str) -> Result<DateTime<Utc>, Report> {
  let parsed = DateTime::<FixedOffset>::parse_from_str(format, date_str)?;
  let utc: DateTime<Utc> = DateTime::from(parsed);
  Ok(utc)
}

pub fn date_from_formats(date_str: &str) -> Result<DateTime<Utc>, Report> {
  const DATE_FORMATS: &[&str] = &[
    "%Y-%m-%d", // 2022-07-23
    "%Y/%m/%d", // 2022/07/23
    "%Y-%m",    // 2022-07
    "%Y/%m",    // 2022/07
  ];
  for format in DATE_FORMATS {
    if let Ok(date) = date_from_format(format, date_str) {
      return Ok(date);
    }
  }
  make_error!("Unknown date format: '{date_str}'")
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

pub fn date_to_year_fraction(date: &DateTime<Utc>) -> f64 {
  let date = date.naive_utc().date();
  let year_start = NaiveDate::from_ymd(date.year(), 1, 1);
  let days_from_year_start = date.signed_duration_since(year_start).num_days();
  let days_in_year = days_in_year(date.year());
  let frac = days_from_year_start as f64 / days_in_year as f64;
  date.year() as f64 + frac
}
