use derive_more::{Display, Error};
use eyre::Report;
use num_traits::ToPrimitive;
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::str::FromStr;
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::make_error;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Display, Serialize, Deserialize)]
#[serde(try_from = "usize", into = "usize")]
pub struct MaxGridPoints(usize);

impl MaxGridPoints {
  pub const MIN: usize = 1_000;

  pub const DEFAULT: usize = 1_000_000;

  pub fn new(points: usize) -> Result<Self, InvalidMaxGridPoints> {
    if points < Self::MIN {
      return Err(InvalidMaxGridPoints { points });
    }
    Ok(Self(points))
  }

  pub fn get(self) -> usize {
    self.0
  }

  #[expect(
    clippy::as_conversions,
    reason = "the limit is far below 2^53, so the conversion to f64 is exact"
  )]
  pub fn point_count(self, required: f64, range: (f64, f64), dx: f64) -> Result<usize, Report> {
    if required.is_nan() {
      return make_error!(
        "Cannot build a grid over [{}, {}] with spacing {}: the number of points is not a number",
        number_text(range.0),
        number_text(range.1),
        number_text(dx)
      );
    }
    if required > self.0 as f64 {
      return Err(Report::new(GridPointLimitExceeded {
        required,
        limit: self,
        range,
        dx,
      }));
    }
    match required.to_usize() {
      Some(points) => Ok(points),
      None => make_error!(
        "Cannot build a grid over [{}, {}] with spacing {}: the number of points is {}",
        number_text(range.0),
        number_text(range.1),
        number_text(dx),
        number_text(required)
      ),
    }
  }
}

impl Default for MaxGridPoints {
  fn default() -> Self {
    Self(Self::DEFAULT)
  }
}

impl TryFrom<usize> for MaxGridPoints {
  type Error = InvalidMaxGridPoints;

  fn try_from(points: usize) -> Result<Self, InvalidMaxGridPoints> {
    Self::new(points)
  }
}

impl From<MaxGridPoints> for usize {
  fn from(limit: MaxGridPoints) -> Self {
    limit.0
  }
}

impl FromStr for MaxGridPoints {
  type Err = InvalidMaxGridPointsText;

  fn from_str(text: &str) -> Result<Self, InvalidMaxGridPointsText> {
    let points = text.trim().parse::<usize>().map_err(|err| InvalidMaxGridPointsText {
      message: format!("'{text}' is not a whole number: {err}"),
    })?;
    Self::new(points).map_err(|err| InvalidMaxGridPointsText {
      message: err.to_string(),
    })
  }
}

impl JsonSchema for MaxGridPoints {
  fn inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("MaxGridPoints")
  }

  fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
    json_schema!({
      "type": "integer",
      "format": "uint",
      "minimum": Self::MIN,
    })
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Display, Error)]
#[display("the grid point limit must be at least {} points, got {points}", MaxGridPoints::MIN)]
pub struct InvalidMaxGridPoints {
  #[error(not(source))]
  pub points: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Display, Error)]
#[display("{message}")]
pub struct InvalidMaxGridPointsText {
  #[error(not(source))]
  message: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Display, Error)]
#[display(
  "A grid over [{}, {}] with spacing {} needs {} the limit of {limit}",
  number_text(range.0),
  number_text(range.1),
  number_text(*dx),
  required_text(*required)
)]
pub struct GridPointLimitExceeded {
  #[error(not(source))]
  pub required: f64,
  #[error(not(source))]
  pub limit: MaxGridPoints,
  #[error(not(source))]
  pub range: (f64, f64),
  #[error(not(source))]
  pub dx: f64,
}

impl GridPointLimitExceeded {
  pub fn required_points(&self) -> Option<usize> {
    required_points(self.required)
  }
}

fn required_points(required: f64) -> Option<usize> {
  required.is_finite().then(|| required.to_usize()).flatten()
}

fn number_text(value: f64) -> String {
  float_to_significant_digits(value, 6)
}

fn required_text(required: f64) -> String {
  required_points(required).map_or_else(
    || "more points than".to_owned(),
    |points| format!("{points} points, more than"),
  )
}
