use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

/// Branch-variance model for root-to-tip statistics.
///
/// Internal branches contribute `variance_factor * branch_length + variance_offset`.
/// Terminal (leaf) branches add a further `variance_offset_leaf`. This mirrors
/// the v0 `TreeRegression` weighting; the defaults `(0, 0, 1)` reproduce v0
/// no-covariation `min_dev` rooting (internal variance 0, leaf variance 1), where
/// the objective is the plain variance of root-to-tip distances.
#[derive(Debug, Clone, Copy, SmartDefault, Serialize, Deserialize)]
pub struct VarianceModel {
  /// Variance contribution proportional to branch length.
  #[default = 0.0]
  pub variance_factor: f64,

  /// Constant variance offset added to every branch.
  #[default = 0.0]
  pub variance_offset: f64,

  /// Additional variance offset added to terminal (leaf) branches.
  #[default = 1.0]
  pub variance_offset_leaf: f64,
}

impl VarianceModel {
  /// Variance of a branch of the given length, excluding the leaf offset.
  ///
  /// Used as the full-branch variance for an edge being split; the split
  /// fractions and the leaf offset are applied by the cost function.
  pub fn branch(&self, branch_length: f64) -> f64 {
    self.variance_factor * branch_length + self.variance_offset
  }

  /// Variance of a terminal branch of the given length, including the leaf offset.
  pub fn leaf_branch(&self, branch_length: f64) -> f64 {
    self.branch(branch_length) + self.variance_offset_leaf
  }
}
