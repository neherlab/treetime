use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MetricsConfig {
  pub(crate) pointwise: PointwiseConfig,
  pub(crate) spatial: SpatialConfig,
  pub(crate) distribution: DistributionConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct PointwiseConfig {
  #[default = 1e-15]
  pub(crate) epsilon: f64,
  #[default = 1e-10]
  pub(crate) log_threshold: f64,
  #[default = 1e-12]
  pub(crate) monotonicity_eta: f64,
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub(crate) abs_tolerances: [f64; 3],
  #[default(_code = "[0.01, 0.001, 0.0001]")]
  pub(crate) rel_tolerances: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct SpatialConfig {
  #[default = 1e-15]
  pub(crate) epsilon: f64,
  #[default = 1e-10]
  pub(crate) log_threshold: f64,
  #[default = 3.0]
  pub(crate) peak_region_radius: f64,
  #[default = 1e-6]
  pub(crate) tail_threshold: f64,
  #[default = 5]
  pub(crate) window_half_width: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct DistributionConfig {
  #[default = 50]
  pub(crate) histogram_bins: usize,
  #[default = 1e-16]
  pub(crate) log_min_value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct ToleranceThresholds {
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub(crate) abs_tolerances: [f64; 3],
  #[default(_code = "[0.01, 0.001, 0.0001]")]
  pub(crate) rel_tolerances: [f64; 3],
  #[default(_code = "[0.999999, 0.9999, 0.99]")]
  pub(crate) r2_thresholds: [f64; 3],
}
