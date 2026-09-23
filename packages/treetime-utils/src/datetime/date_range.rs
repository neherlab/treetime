use crate::datetime::datetime::iso;
use chrono::{DateTime, Duration, NaiveDate, Utc};
use getset::Getters;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, SmartDefault, Clone, Copy, Eq, PartialEq, Serialize, Deserialize, Getters)]
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

  #[allow(
    clippy::unwrap_used,
    clippy::as_conversions,
    reason = "constructs a range from a caller-supplied valid (year, month, day); the year fits i32 and the components form a valid calendar instant"
  )]
  pub fn from_ymd(begin: (u32, u32, u32), end: (u32, u32, u32)) -> Self {
    let begin = NaiveDate::from_ymd_opt(begin.0 as i32, begin.1, begin.2)
      .unwrap()
      .and_hms_opt(0, 0, 0)
      .unwrap()
      .and_utc();

    let end = NaiveDate::from_ymd_opt(end.0 as i32, end.1, end.2)
      .unwrap()
      .and_hms_nano_opt(23, 59, 59, 999_999_999)
      .unwrap()
      .and_utc();

    Self::new(begin, end)
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
