use crate::distribution::distribution_formula::DistributionFormula;
use crate::distribution::distribution_negation::distribution_negation_inplace;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::distribution::y_axis_policy::{NegLog, Plain, YAxisPolicy};
use crate::distribution::{distribution_function::DistributionFunction, distribution_negation::distribution_negation};
use crate::io::dates_csv::DateOrRange;
use approx::ulps_eq;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use treetime_utils::make_error;

pub const TIME_LIMIT: f64 = 1e10;
pub const TIME_EPSILON: f64 = 1e-10;

#[must_use]
#[allow(variant_size_differences)]
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
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

  #[allow(clippy::needless_pass_by_value)]
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

  pub fn likely_time(&self) -> Option<f64> {
    match self {
      Self::Empty => {
        None //
      },
      Self::Point(p) => {
        Some(p.t()) //
      },
      Self::Range(r) => {
        Some(f64::midpoint(r.start(), r.end())) //
      },
      Self::Function(f) => {
        f.likely_time() //
      },
      Self::Formula(f) => {
        Some(f.likely_time()) //
      },
    }
  }

  pub fn t(&self) -> Array1<f64> {
    match self {
      Self::Empty => {
        ndarray::array![] //
      },
      Self::Point(p) => {
        ndarray::array![p.t()] //
      },
      Self::Range(r) => {
        ndarray::array![r.start(), r.end()] //
      },
      Self::Function(f) => {
        f.t().to_owned() //
      },
      Self::Formula(f) => {
        ndarray::array![f.t_min(), f.t_max()] //
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

  pub fn negate(&self) -> Self {
    distribution_negation(self)
  }

  pub fn negate_inplace(&mut self) {
    distribution_negation_inplace(self);
  }

  pub fn time_bounds(&self) -> (f64, f64) {
    let t = self.t();
    debug_assert!(!t.is_empty(), "Cannot extract time bounds from empty distribution");
    (t[0], t[t.len() - 1])
  }

  pub fn eval(&self, t: f64) -> Result<f64, Report> {
    match self {
      Self::Function(f) => f.interp(t),
      Self::Formula(f) => f.eval_single(t),
      Self::Point(p) => {
        if ulps_eq!(t, p.t(), max_ulps = 10) {
          Ok(p.amplitude())
        } else {
          make_error!("Cannot evaluate point distribution outside its support")
        }
      },
      Self::Range(r) => {
        if t >= r.start() && t <= r.end() {
          Ok(r.amplitude())
        } else {
          make_error!("Cannot evaluate range distribution outside its support")
        }
      },
      Self::Empty => make_error!("Cannot evaluate empty distribution"),
    }
  }

  pub fn eval_many(&self, t: &Array1<f64>) -> Result<Array1<f64>, Report> {
    match self {
      Self::Function(f) => f.interp_many(t),
      Self::Formula(f) => f.eval_many(t),
      Self::Point(p) => {
        let results = t
          .iter()
          .map(|&ti| {
            if ulps_eq!(ti, p.t(), max_ulps = 10) {
              Ok(p.amplitude())
            } else {
              make_error!("Cannot evaluate point distribution outside its support")
            }
          })
          .collect::<Result<Vec<f64>, Report>>()?;
        Ok(Array1::from(results))
      },
      Self::Range(r) => {
        let results = t
          .iter()
          .map(|&ti| {
            if ti >= r.start() && ti <= r.end() {
              Ok(r.amplitude())
            } else {
              make_error!("Cannot evaluate range distribution outside its support")
            }
          })
          .collect::<Result<Vec<f64>, Report>>()?;
        Ok(Array1::from(results))
      },
      Self::Empty => make_error!("Cannot evaluate empty distribution"),
    }
  }
}

impl Distribution<Plain> {
  pub fn from_date_or_range(date_or_range: &DateOrRange) -> Self {
    match date_or_range {
      DateOrRange::YearFraction(t) => {
        Self::point(*t, 1.0) //
      },
      DateOrRange::YearFractionRange((start, end)) => {
        Self::range((*start, *end), 1.0) //
      },
    }
  }

  /// Returns the maximum value of the distribution.
  /// - Empty: returns 0.0
  /// - Point: returns amplitude
  /// - Range: returns amplitude
  /// - Function: returns max of y values
  /// - Formula: panics (not supported)
  pub fn max_value(&self) -> f64 {
    match self {
      Distribution::Empty => 0.0,
      Distribution::Point(p) => p.amplitude(),
      Distribution::Range(r) => r.amplitude(),
      Distribution::Function(f) => f.y().max().ok().copied().unwrap_or(0.0),
      Distribution::Formula(_) => {
        panic!("max_value not supported for Formula distributions")
      },
    }
  }

  /// Returns a new distribution with all values multiplied by factor.
  pub fn scale_by(&self, factor: f64) -> Self {
    match self {
      Distribution::Empty => Distribution::Empty,
      Distribution::Point(p) => Distribution::point(p.t(), p.amplitude() * factor),
      Distribution::Range(r) => Distribution::range((r.start(), r.end()), r.amplitude() * factor),
      Distribution::Function(f) => f.scale_y(factor).map_or(Distribution::Empty, Distribution::Function),
      Distribution::Formula(_) => {
        panic!("scale_by not supported for Formula distributions")
      },
    }
  }

  /// Returns a new distribution with max value = 1.0.
  /// Divides all values by max_value().
  /// Returns Empty if max_value() <= 0.
  pub fn normalize(&self) -> Self {
    let max_val = self.max_value();
    if max_val <= 0.0 || !max_val.is_finite() {
      return Distribution::Empty;
    }
    self.scale_by(1.0 / max_val)
  }

  pub fn to_neglog(&self) -> Distribution<NegLog> {
    match self {
      Self::Empty => Distribution::Empty,
      Self::Point(p) => Distribution::point(p.t(), NegLog::from_plain(p.amplitude())),
      Self::Range(r) => Distribution::range((r.start(), r.end()), NegLog::from_plain(r.amplitude())),
      Self::Function(f) => {
        let y_neglog = f.y().mapv(NegLog::from_plain);
        Distribution::Function(DistributionFunction::from_grid_fn(
          treetime_grid::GridFn::from_start_dx_values(f.x_min(), f.dx(), y_neglog)
            .expect("Grid construction should not fail for valid input"),
        ))
      },
      Self::Formula(_) => panic!("Conversion to NegLog not implemented for Formula distributions"),
    }
  }
}

impl Distribution<NegLog> {
  pub fn to_plain(&self) -> Distribution<Plain> {
    match self {
      Self::Empty => Distribution::Empty,
      Self::Point(p) => Distribution::point(p.t(), NegLog::to_plain(p.amplitude())),
      Self::Range(r) => Distribution::range((r.start(), r.end()), NegLog::to_plain(r.amplitude())),
      Self::Function(f) => {
        let y_plain = f.y().mapv(NegLog::to_plain);
        Distribution::Function(DistributionFunction::from_grid_fn(
          treetime_grid::GridFn::from_start_dx_values(f.x_min(), f.dx(), y_plain)
            .expect("Grid construction should not fail for valid input"),
        ))
      },
      Self::Formula(f) => {
        let t_min = f.t_min();
        let t_max = f.t_max();
        let f = f.clone();
        let eval_fn = move |t: f64| -> Result<f64, Report> {
          let y = f.eval_single(t)?;
          Ok(NegLog::to_plain(y))
        };
        Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max))
      },
    }
  }
}

pub type DistributionPlain = Distribution<Plain>;
pub type DistributionNegLog = Distribution<NegLog>;
