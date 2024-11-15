use getset::Getters;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// Represents rectangular signal
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Getters)]
#[getset(get = "pub")]
pub struct DistributionRange<T: Clone + Copy + Debug> {
  x: (T, T),
  y: T,
}
