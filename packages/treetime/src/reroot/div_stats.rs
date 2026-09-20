use crate::reroot::traits::RootStats;
use getset::CopyGetters;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Sub};

#[must_use]
#[derive(Debug, Default, Clone, Copy, PartialEq, Serialize, Deserialize, CopyGetters)]
#[getset(get_copy = "pub")]
pub struct DivStats {
  count: f64,
  d_sum: f64,
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
    let denom = 1.0 / (1.0 + variance * self.count);
    Self {
      count: self.count * denom,
      d_sum: (self.d_sum + branch_length * self.count) * denom,
      dsq_sum: self.dsq_sum + 2.0 * branch_length * self.d_sum + branch_length.powi(2) * self.count
        - variance * (self.d_sum + branch_length * self.count).powi(2) * denom,
    }
  }

  fn score(&self) -> f64 {
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
