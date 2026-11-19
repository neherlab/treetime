use crate::DistributionFunction;
use crate::distribution_core::function::costs_coincide;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::mass_domain::log_mass;
use crate::distribution_ops::negate::distribution_negation;
use eyre::Report;
use ndarray::{Array1, array};
use ndarray_stats::QuantileExt;
use ndarray_stats::errors::MinMaxError;
use std::fmt::Debug;
use strum_macros::Display;
use treetime_grid::{BoundaryBehavior, Side, coordinates_coincide};
use treetime_utils::make_error;

#[must_use]
#[derive(Clone, Debug, Default, PartialEq, Display)]
#[strum(serialize_all = "kebab-case")]
pub enum Distribution {
  #[default]
  Empty,
  Point(DistributionPoint),
  Range(DistributionRange),
  Function(DistributionFunction),
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

    if x.len() == 2 && costs_coincide(y[0], y[1]) {
      return Ok(Self::range((x[0], x[1]), y[1]));
    }

    Ok(Self::Function(DistributionFunction::from_arrays(&x, y)?))
  }

  pub const fn is_point(&self) -> bool {
    matches!(self, Self::Point(_))
  }

  pub fn likely_time(&self) -> Result<Option<f64>, Report> {
    match self {
      Self::Empty => Ok(None),
      Self::Point(p) => Ok(Some(p.t())),
      Self::Range(r) => Ok(Some(f64::midpoint(r.start(), r.end()))),
      Self::Function(f) => f.likely_time(),
    }
  }

  pub fn t(&self) -> Array1<f64> {
    match self {
      Self::Empty => {
        array![]
      },
      Self::Point(p) => {
        array![p.t()]
      },
      Self::Range(r) => {
        array![r.start(), r.end()]
      },
      Self::Function(f) => f.t().to_owned(),
    }
  }

  pub fn y(&self) -> Array1<f64> {
    match self {
      Self::Point(p) => array![p.amplitude()],
      Self::Range(r) => array![r.amplitude(), r.amplitude()],
      Self::Function(f) => f.values(),
      Self::Empty => array![],
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
    }
  }

  pub fn support_bounds(&self) -> Option<(f64, f64)> {
    match self {
      Self::Function(f) => Some((
        hard_limit(f.left_extrap(), f.x_min(), f64::NEG_INFINITY),
        hard_limit(f.right_extrap(), f.x_max(), f64::INFINITY),
      )),
      other => other.time_bounds(),
    }
  }

  pub fn eval(&self, t: f64) -> Result<f64, Report> {
    match self {
      Self::Function(f) => f.interp(t),
      Self::Point(p) => {
        if coordinates_coincide(t, p.t(), t.abs().max(p.t().abs())) {
          Ok(p.amplitude())
        } else {
          Ok(f64::INFINITY)
        }
      },
      Self::Range(r) => {
        if t >= r.start() && t <= r.end() {
          Ok(r.amplitude())
        } else {
          Ok(f64::INFINITY)
        }
      },
      Self::Empty => Ok(f64::INFINITY),
    }
  }

  pub(crate) fn with_left_extrap(self, behavior: BoundaryBehavior) -> Self {
    match self {
      Self::Function(f) => Self::Function(f.with_left_extrap(behavior)),
      other => other,
    }
  }

  pub(crate) fn with_right_extrap(self, behavior: BoundaryBehavior) -> Self {
    match self {
      Self::Function(f) => Self::Function(f.with_right_extrap(behavior)),
      other => other,
    }
  }

  pub(crate) fn fit_soft_tail(self, side: Side, n_fit: usize) -> Result<Self, Report> {
    match self {
      Self::Function(function) => Ok(Self::Function(function.fit_soft_tail(side, n_fit)?)),
      other => Ok(other),
    }
  }

  pub fn scale_time(&self, factor: f64) -> Result<Self, Report> {
    if factor.is_nan() || factor <= 0.0 {
      return make_error!("Time scale factor must be positive, got {factor}");
    }
    match self {
      Self::Empty => Ok(Self::Empty),
      Self::Point(p) => Ok(Self::point(p.t() * factor, p.amplitude())),
      Self::Range(r) => Ok(Self::range((r.start() * factor, r.end() * factor), r.amplitude())),
      Self::Function(f) => Ok(Self::Function(f.scale_time(factor)?)),
    }
  }

  pub fn repeak(&self) -> Result<Self, Report> {
    match self {
      Self::Function(f) => Ok(f.repeak()?.map_or(Self::Empty, Self::Function)),
      Self::Point(p) => Ok(defined_amplitude(p.amplitude(), "point")?.map_or(Self::Empty, |_| self.clone())),
      Self::Range(r) => Ok(defined_amplitude(r.amplitude(), "range")?.map_or(Self::Empty, |_| self.clone())),
      Self::Empty => Ok(Self::Empty),
    }
  }

  pub fn log_mass(&self) -> Result<f64, Report> {
    log_mass(self)
  }

  pub fn normalize(&self) -> Result<Self, Report> {
    match self {
      Distribution::Empty => Ok(Distribution::Empty),
      Distribution::Point(p) => {
        Ok(neglog_peak(p.amplitude(), "point")?.map_or(Distribution::Empty, |_| Distribution::point(p.t(), 0.0)))
      },
      Distribution::Range(r) => Ok(
        neglog_peak(r.amplitude(), "range")?
          .map_or(Distribution::Empty, |_| Distribution::range((r.start(), r.end()), 0.0)),
      ),
      Distribution::Function(f) => neglog_function_normalize(f),
    }
  }
}

fn hard_limit(behavior: BoundaryBehavior, edge: f64, unbounded: f64) -> f64 {
  match behavior {
    BoundaryBehavior::Error | BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_) => edge,
    BoundaryBehavior::HardApproach(law) => law.t_hard,
    BoundaryBehavior::Linear(_) => unbounded,
  }
}

fn neglog_function_normalize(function: &DistributionFunction) -> Result<Distribution, Report> {
  let minimum = match function.y().min() {
    Ok(&minimum) => minimum,
    Err(MinMaxError::EmptyInput) => return Ok(Distribution::Empty),
    Err(MinMaxError::UndefinedOrder) => {
      return make_error!(
        "Cannot normalize a distribution on [{}, {}]: its negative log-likelihood values contain NaN",
        function.x_min(),
        function.x_max()
      );
    },
  };
  let context = format!("function on [{}, {}]", function.x_min(), function.x_max());
  Ok(neglog_peak(minimum, &context)?.map_or(Distribution::Empty, |minimum| {
    Distribution::Function(function.shift_y(-minimum).with_amplitude(0.0))
  }))
}

fn defined_amplitude(amplitude: f64, context: &str) -> Result<Option<f64>, Report> {
  if amplitude.is_nan() {
    return make_error!("Cannot re-peak a distribution {context}: its amplitude is NaN");
  }
  Ok(amplitude.is_finite().then_some(amplitude))
}

fn neglog_peak(minimum: f64, context: &str) -> Result<Option<f64>, Report> {
  if minimum.is_finite() {
    Ok(Some(minimum))
  } else if minimum.is_infinite() && minimum.is_sign_positive() {
    Ok(None)
  } else {
    make_error!("Cannot normalize a distribution {context}: its peak negative log-likelihood is {minimum}")
  }
}
