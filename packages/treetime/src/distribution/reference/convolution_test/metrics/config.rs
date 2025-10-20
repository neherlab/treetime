use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MetricsConfig {
  pub pointwise: PointwiseConfig,
  pub spatial: SpatialConfig,
  pub distribution: DistributionConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct PointwiseConfig {
  #[default = 1e-15]
  pub epsilon: f64,
  #[default = 1e-10]
  pub log_threshold: f64,
  #[default = 1e-12]
  pub monotonicity_eta: f64,
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub abs_tolerances: [f64; 3],
  #[default(_code = "[0.01, 0.001, 0.0001]")]
  pub rel_tolerances: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct SpatialConfig {
  #[default = 1e-15]
  pub epsilon: f64,
  #[default = 1e-10]
  pub log_threshold: f64,
  #[default = 3.0]
  pub peak_region_radius: f64,
  #[default = 1e-6]
  pub tail_threshold: f64,
  #[default = 5]
  pub window_half_width: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct DistributionConfig {
  #[default = 50]
  pub histogram_bins: usize,
  #[default = 1e-16]
  pub log_min_value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
pub struct ToleranceThresholds {
  #[default(_code = "[1e-6, 1e-9, 1e-12]")]
  pub abs_tolerances: [f64; 3],
  #[default(_code = "[0.01, 0.001, 0.0001]")]
  pub rel_tolerances: [f64; 3],
  #[default(_code = "[0.999999, 0.9999, 0.99]")]
  pub r2_thresholds: [f64; 3],
}
