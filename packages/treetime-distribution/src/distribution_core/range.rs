use getset::Getters;

#[derive(Debug, PartialEq, Clone, Getters)]
#[getset(get = "pub")]
pub struct DistributionRange {
  range: (f64, f64),
  ampl: f64,
}

impl DistributionRange {
  pub(crate) fn new(x: (f64, f64), y: f64) -> Self {
    DistributionRange { range: x, ampl: y }
  }

  pub(crate) fn start(&self) -> f64 {
    self.range.0
  }

  pub(crate) fn end(&self) -> f64 {
    self.range.1
  }

  pub(crate) fn amplitude(&self) -> f64 {
    self.ampl
  }
}
