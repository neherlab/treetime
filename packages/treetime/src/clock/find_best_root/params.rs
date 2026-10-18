use crate::clock::clock_set::ClockSet;
use smart_default::SmartDefault;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RerootSpec {
  Method(RerootMethod),
  Tips(Vec<GraphNodeKey>),
}

impl Default for RerootSpec {
  fn default() -> Self {
    Self::Method(RerootMethod::default())
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault)]
pub enum RerootMethod {
  #[default]
  LeastSquares,
  MinDev,
  Oldest,
  ClockFilter,
}

#[derive(Copy, Debug, Clone, PartialEq, SmartDefault)]
pub enum RootObjective {
  #[default]
  EstimatedRate,
  FixedRate(f64),
}

impl RootObjective {
  pub(crate) fn score(self, clock_set: &ClockSet) -> f64 {
    match self {
      Self::EstimatedRate => clock_set.chisq(),
      Self::FixedRate(rate) => clock_set.chisq_fixed_rate(rate),
    }
  }

  pub(crate) fn has_positive_rate(self, clock_set: &ClockSet) -> bool {
    match self {
      Self::EstimatedRate => {
        let det = clock_set.determinant();
        det > 0.0 && clock_set.clock_rate(det) > 0.0
      },
      Self::FixedRate(rate) => rate >= 0.0,
    }
  }
}
