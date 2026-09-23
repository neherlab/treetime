use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type DatesMap = BTreeMap<String, Option<DateConstraint>>;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DateConstraint {
  pub raw: String,
  pub value: DateValue,
}

impl DateConstraint {
  pub fn exact(value: f64) -> Self {
    Self {
      raw: value.to_string(),
      value: DateValue::Exact(DateExact { value }),
    }
  }

  #[inline]
  pub fn mean(&self) -> f64 {
    self.value.mean()
  }

  pub fn is_exact(&self) -> bool {
    self.value.is_exact()
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum DateValue {
  Exact(DateExact),
  Uncertain(DateRange),
  Range(DateRange),
}

impl DateValue {
  #[inline]
  fn mean(&self) -> f64 {
    match self {
      DateValue::Exact(d) => d.value,
      DateValue::Uncertain(r) | DateValue::Range(r) => f64::midpoint(r.start, r.end),
    }
  }

  fn is_exact(&self) -> bool {
    matches!(self, DateValue::Exact(_))
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DateExact {
  pub value: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DateRange {
  pub start: f64,
  pub end: f64,
}

impl DateRange {
  pub fn width(&self) -> f64 {
    self.end - self.start
  }

  pub fn contains(&self, value: f64) -> bool {
    value >= self.start && value <= self.end
  }
}
