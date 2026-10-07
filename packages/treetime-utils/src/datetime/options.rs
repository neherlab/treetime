use chrono::{NaiveDate, NaiveTime};
use smart_default::SmartDefault;

#[derive(Debug, Default, Clone)]
pub struct DateParserOptions {
  pub default_time_of_day: TimeOfDay,
}

#[derive(Debug, SmartDefault, Clone, Copy)]
pub enum TimeOfDay {
  Dawn,

  #[default]
  Noon,

  Dusk,

  Custom(NaiveTime),

  CustomFn(fn(&NaiveDate) -> NaiveTime),
}
