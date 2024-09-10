use crate::utils::datetime::datetime::iso;
use chrono::{DateTime, Duration, TimeZone, Utc};
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

  pub fn from_ymd(begin: (u32, u32, u32), end: (u32, u32, u32)) -> Self {
    let begin = Utc.ymd(begin.0 as i32, begin.1, begin.2).and_hms(0, 0, 0);
    let end = Utc
      .ymd(end.0 as i32, end.1, end.2)
      .and_hms_nano(23, 59, 59, 999_999_999);
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
