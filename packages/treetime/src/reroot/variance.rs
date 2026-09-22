use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, Clone, Copy, SmartDefault, Serialize, Deserialize)]
pub struct VarianceModel {
  #[default = 0.0]
  pub variance_factor: f64,

  #[default = 0.0]
  pub variance_offset: f64,

  #[default = 1.0]
  pub variance_offset_leaf: f64,
}

impl VarianceModel {
  pub(crate) fn branch(&self, branch_length: f64) -> f64 {
    self.variance_factor * branch_length + self.variance_offset
  }

  pub(crate) fn leaf_branch(&self, branch_length: f64) -> f64 {
    self.branch(branch_length) + self.variance_offset_leaf
  }
}
