use getset::Getters;
use std::fmt::Debug;

#[derive(Debug, Eq, PartialEq, Clone, Getters)]
#[getset(get = "pub")]
pub struct DistributionRange<T: Clone + Copy + Debug> {
  range: (T, T),
  ampl: T,
}

impl<T: Clone + Copy + Debug> DistributionRange<T> {
  pub(crate) fn new(x: (T, T), y: T) -> Self {
    DistributionRange { range: x, ampl: y }
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
