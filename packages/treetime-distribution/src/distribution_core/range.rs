use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use deser::{Deserialize, Serialize};
use getset::Getters;
use std::fmt::Debug;

#[derive(Debug, Eq, PartialEq, Clone, Getters, Serialize, Deserialize)]
#[getset(get = "pub")]
pub struct DistributionRange<T: Clone + Copy + Debug, Y: YAxisPolicy = Plain> {
  range: (T, T),
  ampl: T,
  #[deser(skip)]
  #[getset(skip)]
  _policy: PolicyMarker<Y>,
}

impl<T: Clone + Copy + Debug, Y: YAxisPolicy> DistributionRange<T, Y> {
  pub(crate) fn new(x: (T, T), y: T) -> Self {
    DistributionRange {
      range: x,
      ampl: y,
      _policy: PolicyMarker::new(),
    }
  }

  pub(crate) fn start(&self) -> T {
    self.range.0
  }

  pub(crate) fn end(&self) -> T {
    self.range.1
  }

  pub(crate) fn amplitude(&self) -> T {
    self.ampl
  }
}
