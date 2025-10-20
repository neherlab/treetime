use chrono::{NaiveDate, NaiveTime};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Represents a time point within the day
#[derive(Debug, SmartDefault, Clone, Copy, Serialize, Deserialize)]
pub enum TimeOfDay {
  /// 00:00:00.000000
  Dawn,

  /// 12:00:00.000000
  #[default]
  Noon,

  /// 23:59:59.999999
  Dusk,

  /// Custom point in the day
  Custom(NaiveTime),

  /// Custom point in the day
  #[serde(skip)]
  CustomFn(fn(&NaiveDate) -> NaiveTime),
}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct DateParserOptions {
  /// If only the date is give, then which time point to chose within the day
  pub default_time_of_day: TimeOfDay,
}
