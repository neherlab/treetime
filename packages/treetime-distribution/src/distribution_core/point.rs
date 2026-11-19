#[derive(Clone, Debug, PartialEq)]
pub struct DistributionPoint {
  t: f64,
  ampl: f64,
}

impl DistributionPoint {
  pub(crate) fn new(t: f64, ampl: f64) -> Self {
    DistributionPoint { t, ampl }
  }

  pub(crate) fn t(&self) -> f64 {
    self.t
  }

  pub(crate) fn amplitude(&self) -> f64 {
    self.ampl
  }
}
