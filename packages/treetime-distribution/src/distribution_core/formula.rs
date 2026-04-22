use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use eyre::Result;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;

const FORMULA_GRID_SIZE: usize = 200;

/// Distribution that evaluates a formula on-demand.
///
/// Stores a closure that computes values at arbitrary points, avoiding discretization
/// errors from pre-computing on a fixed grid.
///
/// The closure must be thread-safe (Fn) and cloneable via Arc.
#[derive(Serialize, Deserialize)]
pub struct DistributionFormula<Y: YAxisPolicy = Plain> {
  /// Formula that evaluates the distribution at a given time point
  #[serde(skip, default = "default_eval_fn")]
  eval_fn: Arc<dyn Fn(f64) -> Result<f64> + Send + Sync>,
  /// Valid time range [t_min, t_max]
  t_min: f64,
  t_max: f64,
  #[serde(skip)]
  _policy: PolicyMarker<Y>,
}

fn default_eval_fn() -> Arc<dyn Fn(f64) -> Result<f64> + Send + Sync> {
  Arc::new(|_t| Ok(0.0))
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

  pub fn eval_single(&self, t: f64) -> Result<f64> {
    (self.eval_fn)(t)
  }

  pub fn eval_many(&self, t: &Array1<f64>) -> Result<Array1<f64>> {
    let mut result = Array1::zeros(t.len());
    for (i, &ti) in t.iter().enumerate() {
      result[i] = self.eval_single(ti)?;
    }
    Ok(result)
  }

  pub fn t_min(&self) -> f64 {
    self.t_min
  }

  pub fn t_max(&self) -> f64 {
    self.t_max
  }

  pub fn likely_time(&self) -> f64 {
    let midpoint = f64::midpoint(self.t_min, self.t_max);
    let n_points = FORMULA_GRID_SIZE;
    let t = Array1::from_shape_fn(n_points, |i| {
      self.t_min + (self.t_max - self.t_min) * (i as f64 / (n_points - 1) as f64)
    });
    let Ok(values) = self.eval_many(&t) else {
      return midpoint;
    };
    match values.argmax() {
      Ok(idx) => t[idx],
      Err(_) => midpoint,
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
