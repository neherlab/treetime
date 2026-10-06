use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;

/// A date as a decimal year and as a calendar day.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub struct YearDate {
  /// The date as a decimal year, for example `2015.47`.
  pub year: f64,
  /// The calendar day of the date, as `YYYY-MM-DD`.
  pub date: String,
}

impl YearDate {
  pub fn new(year: f64) -> Self {
    Self {
      year,
      date: year_fraction_to_datestring(year),
    }
  }
}
