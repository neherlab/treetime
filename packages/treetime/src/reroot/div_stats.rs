use crate::reroot::traits::RootStats;
use getset::CopyGetters;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Sub};

/// Sufficient statistics for the variance of root-to-tip distances.
///
/// These are the divergence components of the clock regression statistics, with
/// `count` (sum of inverse branch variances over tips) replacing the
/// date-gated `norm`: every tip contributes regardless of whether it has a
/// date. The score is the weighted variance of root-to-tip distances, which
/// equals the v0 `min_dev` objective in the date-free regime (least-squares
/// regression with a fixed zero slope and no date variation). Minimizing it
/// places the root so that tips are as equidistant as possible.
#[must_use]
#[derive(Debug, Default, Clone, Copy, PartialEq, Serialize, Deserialize, CopyGetters)]
#[getset(get_copy = "pub")]
pub struct DivStats {
  /// Sum of inverse branch variances over contributing tips (`Σ 1/σ²`).
  count: f64,
  /// Variance-weighted sum of root-to-tip distances (`Σ d/σ²`).
  d_sum: f64,
  /// Variance-weighted sum of squared root-to-tip distances (`Σ d²/σ²`).
  dsq_sum: f64,
}

impl DivStats {
  pub fn new(count: f64, d_sum: f64, dsq_sum: f64) -> Self {
    Self { count, d_sum, dsq_sum }
  }
}

impl RootStats for DivStats {
  fn leaf(_time: Option<f64>, branch_length: f64, variance: f64) -> Self {
    Self {
      count: 1.0 / variance,
      d_sum: branch_length / variance,
      dsq_sum: branch_length.powi(2) / variance,
    }
  }

  fn propagate(&self, branch_length: f64, variance: f64) -> Self {
    // ClockSet::propagate_averages with all time terms (t_sum, tsq_sum, dt_sum)
    // identically zero. `count` plays the role of `norm`.
    let denom = 1.0 / (1.0 + variance * self.count);
    Self {
      count: self.count * denom,
      d_sum: (self.d_sum + branch_length * self.count) * denom,
      dsq_sum: self.dsq_sum + 2.0 * branch_length * self.d_sum + branch_length.powi(2) * self.count
        - variance * (self.d_sum + branch_length * self.count).powi(2) * denom,
    }
  }

  fn score(&self) -> f64 {
    // Weighted variance of root-to-tip distances: (S_dd·n − S_d²) / (2n).
    // Equals ClockSet::chisq with the time-regression term removed. `count` is
    // always positive (every tip contributes 1/variance > 0), so this is finite.
    (self.dsq_sum * self.count - self.d_sum.powi(2)) / (2.0 * self.count)
  }
}

impl Add for DivStats {
  type Output = Self;

  fn add(self, rhs: Self) -> Self {
    Self {
      count: self.count + rhs.count,
      d_sum: self.d_sum + rhs.d_sum,
      dsq_sum: self.dsq_sum + rhs.dsq_sum,
    }
  }
}

impl Sub for DivStats {
  type Output = Self;

  fn sub(self, rhs: Self) -> Self {
    Self {
      count: self.count - rhs.count,
      d_sum: self.d_sum - rhs.d_sum,
      dsq_sum: self.dsq_sum - rhs.dsq_sum,
    }
  }
}
