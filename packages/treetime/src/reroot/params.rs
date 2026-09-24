use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, Clone, Serialize, Deserialize, SmartDefault)]
#[serde(default)]
pub struct BrentParams {
  #[default = 50]
  pub(crate) brent_max_iters: usize,
  #[default = 1e-6]
  pub(crate) brent_tolerance: f64,
}
