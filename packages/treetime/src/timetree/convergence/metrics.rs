use serde::{Deserialize, Serialize};

/// Tracks convergence metrics across timetree optimization iterations.
///
/// Records log-likelihood components and change counts to monitor convergence:
/// - Sequence changes (n_diff) should approach zero
/// - Polytomies resolved (n_resolved) should stabilize
/// - Log-likelihoods should increase or stabilize
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceMetrics {
  /// Number of ancestral sequence changes in this iteration
  pub n_diff: usize,
  /// Number of polytomies resolved in this iteration
  pub n_resolved: usize,
  /// Sequence log-likelihood (log-probability of sequences given the tree and substitution model)
  pub log_lh_seq: Option<f64>,
  /// Positional log-likelihood (log-probability of node positions on the time axis)
  pub log_lh_pos: Option<f64>,
  /// Coalescent log-likelihood (population genetic log-prior on node times)
  pub log_lh_coal: Option<f64>,
  /// Total log-likelihood (sum of available components; absent components are excluded)
  pub log_lh_total: Option<f64>,
}

impl ConvergenceMetrics {
  pub fn has_converged(&self) -> bool {
    self.n_diff == 0 && self.n_resolved == 0
  }
}
