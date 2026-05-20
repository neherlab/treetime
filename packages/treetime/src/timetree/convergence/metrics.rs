use serde::{Deserialize, Serialize};

/// Tracks convergence metrics across timetree optimization iterations.
///
/// Records likelihood components and change counts to monitor convergence:
/// - Sequence changes (n_diff) should approach zero
/// - Polytomies resolved (n_resolved) should stabilize
/// - Likelihoods should increase or stabilize
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceMetrics {
  /// Number of ancestral sequence changes in this iteration
  pub n_diff: usize,
  /// Number of polytomies resolved in this iteration
  pub n_resolved: usize,
  /// Sequence likelihood (probability of observing sequences given tree and substitution model)
  pub lh_seq: Option<f64>,
  /// Positional likelihood (probability of node positions on time axis)
  pub lh_pos: Option<f64>,
  /// Coalescent likelihood (population genetic prior on node times)
  pub lh_coal: Option<f64>,
  /// Total likelihood (sum of available components; absent components are excluded)
  pub lh_total: Option<f64>,
}

impl ConvergenceMetrics {
  pub fn has_converged(&self) -> bool {
    self.n_diff == 0 && self.n_resolved == 0
  }
}
