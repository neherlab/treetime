use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::Plain;
use serde::{Deserialize, Serialize};

/// Threshold for detecting drift from normalized state (max=1.0).
const NORMALIZATION_DRIFT_THRESHOLD: f64 = 1e-10;

/// A distribution decomposed into scale + normalized shape.
///
/// Representation: P(x) = exp(log_scale) * inner(x)
/// Invariant: max(inner) = 1.0
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ScaledDistribution {
  log_scale: f64,
  inner: Distribution<Plain>,
}

impl ScaledDistribution {
  /// Create from log_scale and already-normalized inner distribution.
  /// Caller must ensure inner is normalized (max = 1.0).
  pub fn from_parts(log_scale: f64, inner: Distribution<Plain>) -> Self {
    Self { log_scale, inner }
  }

  /// Create from a plain distribution by extracting scale.
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

  /// Convert back to plain distribution.
  pub fn to_plain(&self) -> Distribution<Plain> {
    if self.log_scale.is_finite() {
      self.inner.scale_by(self.log_scale.exp())
    } else {
      Distribution::Empty
    }
  }

  /// The log of the scale factor.
  pub fn log_scale(&self) -> f64 {
    self.log_scale
  }

  /// The normalized inner distribution (max = 1.0).
  pub fn inner(&self) -> &Distribution<Plain> {
    &self.inner
  }

  /// The peak value: exp(log_scale).
  pub fn peak_value(&self) -> f64 {
    self.log_scale.exp()
  }

  /// Whether this distribution is empty.
  pub fn is_empty(&self) -> bool {
    matches!(self.inner, Distribution::Empty)
  }

  /// Re-normalize if the invariant has drifted.
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
