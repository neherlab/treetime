use crate::policy::PolicyMarker;
use std::fmt::Debug;
use treetime_grid::{Plain, YAxisPolicy};

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
