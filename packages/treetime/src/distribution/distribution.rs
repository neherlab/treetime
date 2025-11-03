use crate::distribution::distribution_negation::distribution_negation_inplace;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::distribution::{distribution_function::DistributionFunction, distribution_negation::distribution_negation};
use crate::io::dates_csv::DateOrRange;
use approx::ulps_eq;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

#[must_use]
#[allow(variant_size_differences)]
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Distribution {
  #[default]
  Empty,
  Point(DistributionPoint<f64>),
  Range(DistributionRange<f64>),
  Function(DistributionFunction<f64>),
}

impl Distribution {
  pub fn empty() -> Self {
    Self::Empty
  }

  pub fn point(x: f64, y: f64) -> Self {
    Self::Point(DistributionPoint::new(x, y))
  }

  pub fn range((x1, x2): (f64, f64), y: f64) -> Self {
    Self::Range(DistributionRange::new((x1, x2), y))
  }

  pub fn function(x: Array1<f64>, y: Array1<f64>) -> Result<Self, Report> {
    assert_eq!(x.shape(), y.shape());

    if x.is_empty() {
      return Ok(Self::empty());
    }

    if x.len() == 1 {
      return Ok(Self::point(x[0], y[0]));
    }

    if x.len() == 2 && ulps_eq!(y[0], y[1], max_ulps = 10) {
      return Ok(Self::range((x[0], x[1]), y[1]));
    }

    Ok(Self::Function(DistributionFunction::new(x, y)?))
  }

  pub fn from_date_or_range(date_or_range: &DateOrRange) -> Self {
    match date_or_range {
      DateOrRange::YearFraction(t) => Self::point(*t, 1.0),
      DateOrRange::YearFractionRange((start, end)) => Self::range((*start, *end), 1.0),
    }
  }

  pub fn likely_time(&self) -> Option<f64> {
    match self {
      Self::Point(p) => Some(p.t()),
      Self::Range(r) => Some(f64::midpoint(r.start(), r.end())),
      Self::Function(f) => f.likely_time(),
      Self::Empty => None,
    }
  }

  pub fn negate(&self) -> Self {
    distribution_negation(self)
  }

  pub fn negate_inplace(&mut self) {
    distribution_negation_inplace(self);
  }
}
