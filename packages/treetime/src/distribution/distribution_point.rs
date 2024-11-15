use getset::Getters;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Represents a single-point spike or delta-function
/// https://en.wikipedia.org/wiki/Dirac_delta_function
#[derive(Debug, Serialize, Deserialize, Getters)]
#[getset(get = "pub")]
pub struct DistributionPoint<T: Clone + Copy + Debug> {
  x: T,
  y: T,
}
