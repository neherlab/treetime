use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use std::fmt::Debug;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DistributionPoint<T: Clone + Copy + Debug, Y: YAxisPolicy = Plain> {
  t: T,
  ampl: T,
  _policy: PolicyMarker<Y>,
}

impl<T: Clone + Copy + Debug, Y: YAxisPolicy> DistributionPoint<T, Y> {
  pub(crate) fn new(t: T, ampl: T) -> Self {
    DistributionPoint {
      t,
      ampl,
      _policy: PolicyMarker::new(),
    }
  }

  pub(crate) fn t(&self) -> T {
    self.t
  }

  pub(crate) fn amplitude(&self) -> T {
    self.ampl
  }
}
