use std::ops::{Add, Sub};

/// Sufficient statistics for scoring a candidate root position.
///
/// Captures the message-passing pattern shared by all root-scoring objectives:
/// per-tip contributions are accumulated up the tree (`leaf`), pushed across
/// branches (`propagate`), and combined across independent subtrees (`Add`). The
/// complementary "rest of tree" message at an internal node is recovered by
/// subtracting a child's contribution from the node aggregate (`Sub`), mirroring
/// the forward regression pass. `score` is the scalar objective the search
/// minimizes (lower is better).
pub trait RootStats: Clone + Default + Add<Output = Self> + Sub<Output = Self> + Send + Sync {
  /// Contribution of a tip toward its parent across a branch of the given length
  /// and variance. `time` is the tip date when available; objectives that do not
  /// use dates ignore it.
  fn leaf(time: Option<f64>, branch_length: f64, variance: f64) -> Self;

  /// Push accumulated statistics across a branch of the given length and variance,
  /// returning the statistics as seen from the other end of the branch.
  #[must_use]
  fn propagate(&self, branch_length: f64, variance: f64) -> Self;

  /// Scalar objective value for these statistics. The search minimizes this.
  fn score(&self) -> f64;
}
