use crate::Distribution;
use crate::policy::Plain;
use serde::{Deserialize, Serialize};

pub(crate) const NORMALIZATION_DRIFT_THRESHOLD: f64 = 1e-10;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ScaledDistribution {
  log_scale: f64,
  inner: Distribution<Plain>,
}

impl ScaledDistribution {
  pub fn from_parts(log_scale: f64, inner: Distribution<Plain>) -> Self {
    Self { log_scale, inner }
  }

  pub fn from_plain(dist: &Distribution<Plain>) -> Self {
    let max_val = dist.max_value();
    if max_val <= 0.0 || !max_val.is_finite() {
      return Self {
        log_scale: f64::NEG_INFINITY,
        inner: Distribution::Empty,
      };
    }
    Self {
      log_scale: max_val.ln(),
      inner: dist.normalize(),
    }
  }

  pub fn to_plain(&self) -> Distribution<Plain> {
    if self.log_scale.is_finite() {
      self.inner.scale_by(self.log_scale.exp())
    } else {
      Distribution::Empty
    }
  }

  pub fn log_scale(&self) -> f64 {
    self.log_scale
  }

  pub fn inner(&self) -> &Distribution<Plain> {
    &self.inner
  }

  pub fn peak_value(&self) -> f64 {
    self.log_scale.exp()
  }

  pub fn is_empty(&self) -> bool {
    matches!(self.inner, Distribution::Empty)
  }

  pub fn renormalize(&mut self) {
    let max_val = self.inner.max_value();
    if max_val <= 0.0 || !max_val.is_finite() {
      self.log_scale = f64::NEG_INFINITY;
      self.inner = Distribution::Empty;
      return;
    }
    if (max_val - 1.0).abs() > NORMALIZATION_DRIFT_THRESHOLD {
      self.log_scale += max_val.ln();
      self.inner = self.inner.normalize();
    }
  }
}

impl Default for ScaledDistribution {
  fn default() -> Self {
    Self {
      log_scale: f64::NEG_INFINITY,
      inner: Distribution::Empty,
    }
  }
}
