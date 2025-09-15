use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use approx::ulps_eq;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

#[allow(variant_size_differences)]
#[derive(Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Distribution {
  #[default]
  Empty,
  Point(DistributionPoint<f64>),
  Range(DistributionRange<f64>),

  #[serde(skip)]
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
}
