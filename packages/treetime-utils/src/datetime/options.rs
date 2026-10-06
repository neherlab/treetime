use chrono::{NaiveDate, NaiveTime};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, Default, Clone, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct DateParserOptions {
  pub default_time_of_day: TimeOfDay,
}

#[derive(Debug, SmartDefault, Clone, Copy, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum TimeOfDay {
  Dawn,

  #[default]
  Noon,

  Dusk,

  Custom(NaiveTime),

  #[serde(skip)]
  #[deser(skip)]
  CustomFn(fn(&NaiveDate) -> NaiveTime),
}
