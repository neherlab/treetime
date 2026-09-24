use crate::DistributionFunction;
use crate::distribution_core::formula::DistributionFormula;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::negate::distribution_negation;
use crate::policy::{NegLog, Plain, YAxisPolicy};
use approx::ulps_eq;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use strum_macros::Display;
use treetime_grid::{BoundaryBehavior, Side};

pub(crate) const TIME_LIMIT: f64 = 1e10;
const FORMULA_GRID_SIZE: usize = 200;

#[must_use]
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, Display)]
#[strum(serialize_all = "kebab-case")]
#[serde(rename_all = "kebab-case")]
pub enum Distribution<Y: YAxisPolicy = Plain> {
  #[default]
  Empty,
  Point(DistributionPoint<f64, Y>),
  Range(DistributionRange<f64, Y>),
  Function(DistributionFunction<f64, Y>),
  Formula(DistributionFormula<Y>),
}

impl<Y: YAxisPolicy> Distribution<Y> {
  pub fn empty() -> Self {
    Self::Empty
  }

  pub fn point(x: f64, y: f64) -> Self {
    Self::Point(DistributionPoint::new(x, y))
  }

  pub fn range((x1, x2): (f64, f64), y: f64) -> Self {
    Self::Range(DistributionRange::new((x1, x2), y))
  }

  #[expect(
    clippy::needless_pass_by_value,
    reason = "the constructor stores the arrays in the distribution it returns"
  )]
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

    Ok(Self::Function(DistributionFunction::from_arrays(&x, y)?))
  }

  pub fn constant(amplitude: f64) -> Self {
    Distribution::range((-TIME_LIMIT, TIME_LIMIT), amplitude)
  }

  pub const fn is_point(&self) -> bool {
    matches!(self, Self::Point(_))
  }

  pub fn likely_time(&self) -> Option<f64> {
    match self {
      Self::Empty => None,
      Self::Point(p) => Some(p.t()),
      Self::Range(r) => Some(f64::midpoint(r.start(), r.end())),
      Self::Function(f) => f.likely_time(),
      Self::Formula(f) => Some(f.likely_time()),
    }
  }

  pub fn t(&self) -> Array1<f64> {
    match self {
      Self::Empty => {
        ndarray::array![]
      },
      Self::Point(p) => {
        ndarray::array![p.t()]
      },
      Self::Range(r) => {
        ndarray::array![r.start(), r.end()]
      },
      Self::Function(f) => f.t().to_owned(),
      Self::Formula(f) => {
        ndarray::array![f.t_min(), f.t_max()]
      },
    }
  }

  pub fn y(&self) -> Array1<f64> {
    match self {
      Self::Point(p) => ndarray::array![p.amplitude()],
      Self::Range(r) => ndarray::array![r.amplitude(), r.amplitude()],
      Self::Function(f) => f.y().clone(),
      Self::Formula(f) => {
        let t = ndarray::array![f.t_min(), f.t_max()];
        f.eval_many(&t).unwrap_or_else(|_| ndarray::array![0.0, 0.0])
      },
      Self::Empty => ndarray::array![],
    }
  }

  pub fn negate(&self) -> Result<Self, Report> {
    distribution_negation(self)
  }

  pub fn time_bounds(&self) -> Option<(f64, f64)> {
    match self {
      Self::Empty => None,
      Self::Point(p) => Some((p.t(), p.t())),
      Self::Range(r) => Some((r.start(), r.end())),
      Self::Function(f) => Some((f.x_min(), f.x_max())),
      Self::Formula(f) => Some((f.t_min(), f.t_max())),
    }
  }

  pub fn eval(&self, t: f64) -> Result<f64, Report> {
    match self {
      Self::Function(f) => f.interp(t),
      Self::Formula(f) => f.eval_single(t),
      Self::Point(p) => {
        if ulps_eq!(t, p.t(), max_ulps = 10) {
          Ok(p.amplitude())
        } else {
          Ok(Y::probability_zero())
        }
      },
      Self::Range(r) => {
        if t >= r.start() && t <= r.end() {
          Ok(r.amplitude())
        } else {
          Ok(Y::probability_zero())
        }
      },
      Self::Empty => Ok(Y::probability_zero()),
    }
  }

  pub(crate) fn with_left_extrap(self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    match self {
      Self::Function(f) => Ok(Self::Function(f.with_left_extrap(behavior)?)),
      other => Ok(other),
    }
  }

  pub(crate) fn with_right_extrap(self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    match self {
      Self::Function(f) => Ok(Self::Function(f.with_right_extrap(behavior)?)),
      other => Ok(other),
    }
  }

  pub(crate) fn fit_soft_tail(self, side: Side, n_fit: usize) -> Result<Self, Report> {
    match self {
      Self::Function(function) => Ok(Self::Function(function.fit_soft_tail(side, n_fit)?)),
      other => Ok(other),
    }
  }
}

impl Distribution<Plain> {}

impl Distribution<NegLog> {
  pub fn normalize(&self) -> Self {
    match self {
      Distribution::Empty => Distribution::Empty,
      Distribution::Point(p) => {
        if p.amplitude().is_finite() {
          Distribution::point(p.t(), 0.0)
        } else {
          Distribution::Empty
        }
      },
      Distribution::Range(r) => {
        if r.amplitude().is_finite() {
          Distribution::range((r.start(), r.end()), 0.0)
        } else {
          Distribution::Empty
        }
      },
      Distribution::Function(f) => neglog_function_normalize(f),
      Distribution::Formula(f) => {
        discretize_formula(f).map_or(Distribution::Empty, |df| neglog_function_normalize(&df))
      },
    }
  }
}

fn neglog_function_normalize(function: &DistributionFunction<f64, NegLog>) -> Distribution<NegLog> {
  let Some(minimum) = function.y().min().ok().copied().filter(|minimum| minimum.is_finite()) else {
    return Distribution::Empty;
  };
  Distribution::Function(function.shift_y(-minimum))
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn discretize_formula<Y: YAxisPolicy>(f: &DistributionFormula<Y>) -> Result<DistributionFunction<f64, Y>, Report> {
  let n_points = FORMULA_GRID_SIZE;
  let t = Array1::from_shape_fn(n_points, |i| {
    f.t_min() + (f.t_max() - f.t_min()) * (i as f64 / (n_points - 1) as f64)
  });
  let values = f.eval_many(&t)?;
  DistributionFunction::from_range_values((f.t_min(), f.t_max()), values)
}
