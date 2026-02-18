use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use eyre::Result;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

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
    unimplemented!(
      "likely_time() not available for DistributionFormula: \
       finding the peak requires discretization. \
       Use discretization (DistributionFunction) if you need likely_time()"
    )
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

impl<Y: YAxisPolicy> std::fmt::Debug for DistributionFormula<Y> {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    f.write_str("DistributionFormula")
  }
}

impl<Y: YAxisPolicy> PartialEq for DistributionFormula<Y> {
  fn eq(&self, other: &Self) -> bool {
    self.t_min == other.t_min && self.t_max == other.t_max
  }
}
