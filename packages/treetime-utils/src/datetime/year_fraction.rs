use crate::cast::{ToF64Exact, f64_to_i32_saturating};
use crate::datetime::date_range::DateRange;
use crate::make_report;
use chrono::{DateTime, Datelike, NaiveDate, NaiveTime, TimeDelta, TimeZone, Utc};
use eyre::Report;
use std::ops::RangeInclusive;
use std::time::Duration;

const SECONDS_PER_DAY: u32 = 24 * 60 * 60;

const NANOSECONDS_PER_SECOND: f64 = 1e9;

const FOUR_DIGIT_YEARS: RangeInclusive<i32> = 0..=9999;

pub fn date_range_to_year_fraction_range(date_range: &DateRange) -> (f64, f64) {
  let begin = date_to_year_fraction(date_range.begin());
  let end = date_to_year_fraction(date_range.end());
  (begin, end)
}

pub fn date_to_year_fraction(date: &DateTime<Utc>) -> f64 {
  let frac = (date.ordinal().to_f64_exact() - 0.5) / days_in_year(date.date_naive()).to_f64_exact();
  date.year().to_f64_exact() + frac
}

pub fn year_fraction_at_month_start(year: i32, month: u32) -> Result<f64, Report> {
  let first_day = NaiveDate::from_ymd_opt(year, month, 1)
    .ok_or_else(|| make_report!("Month {month} of year {year} is not a calendar month"))?;
  Ok(year.to_f64_exact() + first_day.ordinal0().to_f64_exact() / days_in_year(first_day).to_f64_exact())
}

pub fn year_fraction_to_datestring(year_fraction: f64) -> Option<String> {
  year_fraction_to_date(year_fraction).map(|date| calendar_date_string(&date))
}

pub fn year_fraction_days_between(from: f64, to: f64) -> Option<f64> {
  let span = year_fraction_to_date(to)? - year_fraction_to_date(from)?;
  let seconds = span.num_seconds().to_f64_exact() + f64::from(span.subsec_nanos()) / NANOSECONDS_PER_SECOND;
  Some(seconds / f64::from(SECONDS_PER_DAY))
}

pub fn year_fraction_to_date(year_fraction: f64) -> Option<DateTime<Utc>> {
  if !year_fraction.is_finite() {
    return None;
  }
  let year = year_fraction.floor();
  let first_day = NaiveDate::from_ymd_opt(f64_to_i32_saturating(year), 1, 1)?;
  let seconds_in_year = (days_in_year(first_day) * SECONDS_PER_DAY).to_f64_exact();
  let since_year_start = Duration::try_from_secs_f64((year_fraction - year) * seconds_in_year).ok()?;
  let since_year_start = TimeDelta::from_std(since_year_start).ok()?;
  Utc
    .from_utc_datetime(&first_day.and_time(NaiveTime::MIN))
    .checked_add_signed(since_year_start)
}

fn calendar_date_string(date: &DateTime<Utc>) -> String {
  let year = date.year();
  if FOUR_DIGIT_YEARS.contains(&year) {
    date.format("%Y-%m-%d").to_string()
  } else {
    format!("{year:+06}-{}", date.format("%m-%d"))
  }
}

fn days_in_year(date: NaiveDate) -> u32 {
  if date.leap_year() { 366 } else { 365 }
}
