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

/// Optimization method selection
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum OptimizationMethod {
  /// Grid search with equally-spaced evaluation points
  #[default]
  Grid,
  /// Brent's method for robust 1D optimization
  Brent,
  /// Golden section search optimization
  #[cfg_attr(feature = "clap", clap(name = "golden-section"))]
  GoldenSection,
}

impl BranchPointOptimizationParams {
  /// Create Brent's method with default parameters
  pub fn brent() -> Self {
    Self::Brent(BrentParams::default())
  }

  /// Create Brent's method with custom parameters
  pub fn brent_with(params: BrentParams) -> Self {
    Self::Brent(params)
  }

  /// Create golden section search with default parameters
  pub fn golden_section() -> Self {
    Self::GoldenSection(GoldenSectionParams::default())
  }

  /// Create golden section search with custom parameters
  pub fn golden_section_with(params: GoldenSectionParams) -> Self {
    Self::GoldenSection(params)
  }

  /// Create grid search with default parameters
  pub fn grid() -> Self {
    Self::Grid(GridSearchParams::default())
  }

  /// Create grid search with custom parameters
  pub fn grid_with(params: GridSearchParams) -> Self {
    Self::Grid(params)
  }
}

/// Configuration for grid search optimization
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GridSearchParams {
  /// Number of equally-spaced points to evaluate (grid method only)
  #[cfg_attr(feature = "clap", clap(long = "branch-split-grid-n-points", default_value_t = GridSearchParams::default().n_points))]
  #[default = 11]
  pub n_points: usize,
}

/// Configuration for Brent's method optimization
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct BrentParams {
  /// Maximum number of iterations for Brent's method
  #[cfg_attr(feature = "clap", clap(long = "branch-split-brent-max-iters", default_value_t = BrentParams::default().brent_max_iters))]
  #[default = 50]
  pub brent_max_iters: usize,
  /// Convergence tolerance for Brent's method
  #[cfg_attr(feature = "clap", clap(long = "branch-split-brent-tolerance", default_value_t = BrentParams::default().brent_tolerance))]
  #[default = 1e-12]
  pub brent_tolerance: f64,
}

/// Configuration for golden section search optimization
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GoldenSectionParams {
  /// Maximum number of iterations for golden section search
  #[cfg_attr(feature = "clap", clap(long = "branch-split-golden-max-iters", default_value_t = GoldenSectionParams::default().golden_max_iters))]
  #[default = 50]
  pub golden_max_iters: usize,
  /// Convergence tolerance for golden section search
  #[cfg_attr(feature = "clap", clap(long = "branch-split-golden-tolerance", default_value_t = GoldenSectionParams::default().golden_tolerance))]
  #[default = 1e-12]
  pub golden_tolerance: f64,
}
