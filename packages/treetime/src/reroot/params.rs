use crate::clock::clock_set::ClockSet;
use crate::reroot::orchestrate::RerootTopologyParams;
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

#[derive(Clone, Debug, SmartDefault)]
pub struct RerootParams {
  pub(crate) spec: RerootSpec,

  pub(crate) objective: RootObjective,

  pub(crate) topology: RerootTopologyParams,

  #[default = true]
  pub(crate) force_positive_rate: bool,
}

impl RerootParams {
  #[must_use]
  pub(crate) fn new(spec: RerootSpec, force_positive_rate: bool) -> Self {
    Self {
      spec,
      force_positive_rate,
      ..Self::default()
    }
  }

  #[must_use]
  pub(crate) fn with_clock_rate(self, clock_rate: Option<f64>) -> Self {
    match clock_rate {
      Some(rate) => Self {
        objective: RootObjective::FixedRate(rate),
        ..self
      },
      None => self,
    }
  }
}
