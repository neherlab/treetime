use crate::distribution::y_axis_policy::{Plain, PolicyMarker, YAxisPolicy};
use getset::Getters;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Rectangular signal
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize, Getters)]
#[getset(get = "pub")]
pub struct DistributionRange<T: Clone + Copy + Debug, Y: YAxisPolicy = Plain> {
  range: (T, T),
  ampl: T,
  #[serde(skip)]
  #[getset(skip)]
  _policy: PolicyMarker<Y>,
}

impl<T: Clone + Copy + Debug, Y: YAxisPolicy> DistributionRange<T, Y> {
  pub fn new(x: (T, T), y: T) -> Self {
    DistributionRange {
      range: x,
      ampl: y,
      _policy: PolicyMarker::new(),
    }
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
