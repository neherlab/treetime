use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Configuration for Brent's-method split optimization.
#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[serde(default)]
pub struct BrentParams {
  #[default = 50]
  pub brent_max_iters: usize,
  #[default = 1e-6]
  pub brent_tolerance: f64,
}
