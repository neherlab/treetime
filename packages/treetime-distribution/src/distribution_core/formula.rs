use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use eyre::{Result, WrapErr};
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use ndarray_stats::errors::MinMaxError;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;
use treetime_utils::{make_error, make_internal_error};

const FORMULA_GRID_SIZE: usize = 200;

#[derive(Serialize, Deserialize)]
pub struct DistributionFormula<Y: YAxisPolicy = Plain> {
  #[serde(skip, default = "default_eval_fn")]
  eval_fn: Arc<dyn Fn(f64) -> Result<f64> + Send + Sync>,
  t_min: f64,
  t_max: f64,
  #[serde(skip)]
  _policy: PolicyMarker<Y>,
}

impl<Y: YAxisPolicy> DistributionFormula<Y> {
  pub fn new<F>(eval_fn: F, t_min: f64, t_max: f64) -> Self
  where
    F: Fn(f64) -> Result<f64> + Send + Sync + 'static,
  {
    Self {
      eval_fn: Arc::new(eval_fn),
      t_min,
      t_max,
      _policy: PolicyMarker::new(),
    }
  }

  pub(crate) fn eval_single(&self, t: f64) -> Result<f64> {
    (self.eval_fn)(t)
  }

  pub(crate) fn eval_many(&self, t: &Array1<f64>) -> Result<Array1<f64>> {
    let mut result = Array1::zeros(t.len());
    for (i, &ti) in t.iter().enumerate() {
      result[i] = self.eval_single(ti)?;
    }
    Ok(result)
  }

  pub(crate) fn t_min(&self) -> f64 {
    self.t_min
  }

  pub(crate) fn t_max(&self) -> f64 {
    self.t_max
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  pub(crate) fn likely_time(&self) -> Result<f64> {
    let n_points = FORMULA_GRID_SIZE;
    let t = Array1::from_shape_fn(n_points, |i| {
      self.t_min + (self.t_max - self.t_min) * (i as f64 / (n_points - 1) as f64)
    });
    let values = self.eval_many(&t).wrap_err_with(|| {
      format!(
        "When finding the most likely time of a formula distribution on [{}, {}]",
        self.t_min, self.t_max
      )
    })?;
    let extremum = if Y::likely_is_maximum() {
      values.argmax()
    } else {
      values.argmin()
    };
    match extremum {
      Ok(idx) => Ok(t[idx]),
      Err(MinMaxError::UndefinedOrder) => make_error!(
        "Cannot find the most likely time of a formula distribution on [{}, {}]: its values contain NaN",
        self.t_min,
        self.t_max
      ),
      Err(MinMaxError::EmptyInput) => make_internal_error!("The formula evaluation grid has no points"),
    }
  }
}

impl<Y: YAxisPolicy> Clone for DistributionFormula<Y> {
  fn clone(&self) -> Self {
    Self {
      eval_fn: Arc::clone(&self.eval_fn),
      t_min: self.t_min,
      t_max: self.t_max,
      _policy: PolicyMarker::new(),
    }
  }
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    handwritten_fmt_impl,
    reason = "the formula holds a closure, which has no Debug representation"
  )
)]
impl<Y: YAxisPolicy> fmt::Debug for DistributionFormula<Y> {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    f.write_str("DistributionFormula")
  }
}

impl<Y: YAxisPolicy> PartialEq for DistributionFormula<Y> {
  fn eq(&self, other: &Self) -> bool {
    self.t_min == other.t_min && self.t_max == other.t_max
  }
}

fn default_eval_fn() -> Arc<dyn Fn(f64) -> Result<f64> + Send + Sync> {
  Arc::new(|_t| Ok(0.0))
}
