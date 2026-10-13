use crate::datetime::datetime::{LAST_NANOSECOND, iso};
use crate::make_report;
use chrono::{DateTime, Duration, NaiveDate, NaiveTime, Utc};
use eyre::{Report, WrapErr};
use getset::Getters;
use smart_default::SmartDefault;

#[derive(Debug, SmartDefault, Clone, Copy, Eq, PartialEq, Getters)]
#[getset(get = "pub")]
pub struct DateRange {
  begin: DateTime<Utc>,
  end: DateTime<Utc>,
}

impl DateRange {
  pub fn new(begin: impl Into<DateTime<Utc>>, end: impl Into<DateTime<Utc>>) -> Self {
    Self {
      begin: begin.into(),
      end: end.into(),
    }
  }

  pub fn from_iso(begin: impl AsRef<str>, end: impl AsRef<str>) -> Self {
    Self::new(iso(begin.as_ref()), iso(end.as_ref()))
  }

  pub fn from_ymd(begin: (u32, u32, u32), end: (u32, u32, u32)) -> Result<Self, Report> {
    let begin = naive_date(begin)?.and_time(NaiveTime::MIN).and_utc();

    let end = naive_date(end)?.and_time(LAST_NANOSECOND).and_utc();

    Ok(Self::new(begin, end))
  }

  pub fn is_valid(&self) -> bool {
    self.end > self.begin
  }

  pub fn len(&self) -> Duration {
    self.end - self.begin
  }

  pub fn mean(&self) -> DateTime<Utc> {
    self.begin + (self.len() / 2)
  }
}

fn naive_date((year, month, day): (u32, u32, u32)) -> Result<NaiveDate, Report> {
  let signed_year = i32::try_from(year).wrap_err_with(|| format!("Year {year} is outside the supported date range"))?;
  NaiveDate::from_ymd_opt(signed_year, month, day)
    .ok_or_else(|| make_report!("Date {year:04}-{month:02}-{day:02} is not a valid calendar date"))
}
