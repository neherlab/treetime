use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Configuration for branch point optimization methods
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub enum BranchPointOptimizationParams {
  /// Grid search with specified number of equally-spaced points
  #[default]
  Grid(GridSearchParams),

  /// Brent's method for robust 1D optimization
  Brent(BrentParams),

  /// Golden section search
  GoldenSection(GoldenSectionParams),
}

impl BranchPointOptimizationParams {
  /// Create Brent's method with default parameters
  pub fn brent() -> Self {
    Self::Brent(BrentParams::default())
  }

  /// Create golden section search with default parameters
  pub fn golden_section() -> Self {
    Self::GoldenSection(GoldenSectionParams::default())
  }

  /// Create grid search with specified number of points
  pub fn grid() -> Self {
    Self::Grid(GridSearchParams::default())
  }
}

/// Configuration for grid search optimization
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct GridSearchParams {
  /// Number of equally-spaced points to evaluate
  #[default = 11]
  pub n_points: usize,
}

/// Configuration for Brent's method optimization
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct BrentParams {
  /// Maximum number of iterations
  #[default = 50]
  pub max_iters: usize,
  /// Convergence tolerance
  #[default = 1e-12]
  pub tolerance: f64,
}

/// Configuration for golden section search optimization
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct GoldenSectionParams {
  /// Maximum number of iterations
  #[default = 50]
  pub max_iters: usize,
  /// Convergence tolerance
  #[default = 1e-12]
  pub tolerance: f64,
}
