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

#[derive(Debug, Clone, SmartDefault)]
pub enum BranchPointOptimizationParams {
  #[default]
  Grid(GridSearchParams),

  Brent(BrentParams),

  GoldenSection(GoldenSectionParams),
}

impl BranchPointOptimizationParams {
  pub fn brent_with(params: BrentParams) -> Self {
    Self::Brent(params)
  }

  pub fn golden_section_with(params: GoldenSectionParams) -> Self {
    Self::GoldenSection(params)
  }

  pub fn grid_with(params: GridSearchParams) -> Self {
    Self::Grid(params)
  }
}

#[derive(Debug, Clone, SmartDefault)]
pub enum OptimizationMethod {
  #[default]
  Grid,
  Brent,
  GoldenSection,
}

#[derive(Debug, Clone, SmartDefault)]
pub struct GridSearchParams {
  #[default = 11]
  pub n_points: usize,
}

#[derive(Debug, Clone, SmartDefault)]
pub struct BrentParams {
  #[default = 50]
  pub brent_max_iters: usize,
  #[default = 1e-12]
  pub brent_tolerance: f64,
}

#[derive(Debug, Clone, SmartDefault)]
pub struct GoldenSectionParams {
  #[default = 50]
  pub golden_max_iters: usize,
  #[default = 1e-12]
  pub golden_tolerance: f64,
}
