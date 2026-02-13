use crate::y_axis_policy::{Plain, PolicyMarker, YAxisPolicy};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Represents a single-point spike or delta-function
/// https://en.wikipedia.org/wiki/Dirac_delta_function
#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct DistributionPoint<T: Clone + Copy + Debug, Y: YAxisPolicy = Plain> {
  t: T,
  ampl: T,
  #[serde(skip)]
  _policy: PolicyMarker<Y>,
}

impl<T: Clone + Copy + Debug, Y: YAxisPolicy> DistributionPoint<T, Y> {
  pub fn new(t: T, ampl: T) -> Self {
    DistributionPoint {
      t,
      ampl,
      _policy: PolicyMarker::new(),
    }
  }

  pub fn t(&self) -> T {
    self.t
  }

  pub fn amplitude(&self) -> T {
    self.ampl
  }
}
