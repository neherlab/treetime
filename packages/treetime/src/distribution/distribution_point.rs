use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Represents a single-point spike or delta-function
/// https://en.wikipedia.org/wiki/Dirac_delta_function
#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct DistributionPoint<T: Clone + Copy + Debug> {
  t: T,
  ampl: T,
}

impl<T: Clone + Copy + Debug> DistributionPoint<T> {
  pub fn new(x: T, y: T) -> Self {
    DistributionPoint { t: x, ampl: y }
  }

  pub fn t(&self) -> T {
    self.t
  }

  pub fn amplitude(&self) -> T {
    self.ampl
  }
}
