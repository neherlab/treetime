use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;

/// A date as a decimal year and as a calendar day.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct YearDate {
  /// The date as a decimal year, for example `2015.47`.
  pub year: f64,
  /// The calendar day of the date: `YYYY-MM-DD` for the years 0 to 9999, and the signed five-digit year of FITS and ISO 8601 for other years, for example `-01910-05-30`. Year 0 is 1 BCE. Absent when the decimal year lies outside the years -262143 to 262142.
  pub date: Option<String>,
}

impl YearDate {
  pub fn new(year: f64) -> Self {
    Self {
      year,
      date: year_fraction_to_datestring(year),
    }
  }
}
