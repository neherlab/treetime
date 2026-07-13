use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Selection of the 1D method used to optimize the split position along a candidate edge.
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub enum BranchPointOptimizationParams {
  /// Grid search with equally-spaced evaluation points.
  #[default]
  Grid(GridSearchParams),

  /// Brent's method.
  Brent(BrentParams),

  /// Golden section search.
  GoldenSection(GoldenSectionParams),
}

impl BranchPointOptimizationParams {
  pub fn grid() -> Self {
    Self::Grid(GridSearchParams::default())
  }

  pub fn brent() -> Self {
    Self::Brent(BrentParams::default())
  }

  pub fn golden_section() -> Self {
    Self::GoldenSection(GoldenSectionParams::default())
  }
}

/// Configuration for grid-search split optimization.
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[serde(default)]
pub struct GridSearchParams {
  /// Number of equally-spaced points evaluated along the edge.
  #[default = 11]
  pub n_points: usize,
}

/// Configuration for Brent's-method split optimization.
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[serde(default)]
pub struct BrentParams {
  #[default = 50]
  pub brent_max_iters: usize,
  #[default = 1e-12]
  pub brent_tolerance: f64,
}

/// Configuration for golden-section split optimization.
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[serde(default)]
pub struct GoldenSectionParams {
  #[default = 50]
  pub golden_max_iters: usize,
  #[default = 1e-12]
  pub golden_tolerance: f64,
}
