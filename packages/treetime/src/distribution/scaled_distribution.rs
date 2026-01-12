use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::Plain;
use serde::{Deserialize, Serialize};

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
    if (max_val - 1.0).abs() > 1e-10 {
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

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_relative_eq;
  use ndarray::array;

  #[test]
  fn test_scaled_distribution_from_plain_round_trip() {
    let original = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]).unwrap();

    let scaled = ScaledDistribution::from_plain(&original);
    let recovered = scaled.to_plain();

    if let (Distribution::Function(orig), Distribution::Function(rec)) = (&original, &recovered) {
      for i in 0..3 {
        assert_relative_eq!(orig.y()[i], rec.y()[i], epsilon = 1e-10);
      }
    } else {
      panic!("Expected Function distributions");
    }
  }

  #[test]
  fn test_scaled_distribution_log_scale() {
    let dist = Distribution::<Plain>::point(1.0, 100.0);
    let scaled = ScaledDistribution::from_plain(&dist);

    assert_relative_eq!(scaled.log_scale(), 100.0_f64.ln());
    assert_relative_eq!(scaled.peak_value(), 100.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_inner_is_normalized() {
    let dist = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]).unwrap();
    let scaled = ScaledDistribution::from_plain(&dist);

    assert_relative_eq!(scaled.inner().max_value(), 1.0);
  }

  #[test]
  fn test_scaled_distribution_empty() {
    let empty = Distribution::<Plain>::Empty;
    let scaled = ScaledDistribution::from_plain(&empty);

    assert!(scaled.is_empty());
    assert!(scaled.log_scale().is_infinite());
  }

  #[test]
  fn test_scaled_distribution_default() {
    let default = ScaledDistribution::default();
    assert!(default.is_empty());
    assert!(default.log_scale().is_infinite());
  }

  #[test]
  fn test_scaled_distribution_renormalize() {
    let dist = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![10.0, 40.0, 20.0]).unwrap();
    let mut scaled = ScaledDistribution::from_plain(&dist);

    let original_log_scale = scaled.log_scale();
    scaled.renormalize();

    assert_relative_eq!(scaled.log_scale(), original_log_scale, epsilon = 1e-10);
    assert_relative_eq!(scaled.inner().max_value(), 1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_scaled_distribution_from_parts() {
    let inner = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![0.25, 1.0, 0.5]).unwrap();
    let log_scale = 2.0_f64.ln();
    let scaled = ScaledDistribution::from_parts(log_scale, inner);

    assert_relative_eq!(scaled.log_scale(), log_scale);
    assert_relative_eq!(scaled.peak_value(), 2.0, epsilon = 1e-10);
  }
}
