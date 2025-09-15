use getset::Getters;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Represents rectangular signal
#[derive(Debug, Eq, PartialEq, Clone, Copy, Serialize, Deserialize, Getters)]
#[getset(get = "pub")]
pub struct DistributionRange<T: Clone + Copy + Debug> {
  range: (T, T),
  ampl: T,
}

impl<T: Clone + Copy + Debug> DistributionRange<T> {
  pub fn new(x: (T, T), y: T) -> Self {
    DistributionRange { range: x, ampl: y }
  }

  pub fn start(&self) -> T {
    self.range.0
  }

  pub fn end(&self) -> T {
    self.range.1
  }

  pub fn amplitude(&self) -> T {
    self.ampl
  }
}
