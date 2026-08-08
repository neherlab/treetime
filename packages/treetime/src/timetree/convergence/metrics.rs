use serde::{Deserialize, Serialize};
use treetime_primitives::LogLh;

/// Largest node-time movement, in years, that still counts as converged.
///
/// Set at the resolution the data carry rather than at the resolution the solver can report:
/// sampling dates are rarely recorded finer than a day, and the biological replication cycle is
/// itself a few days, so movement below ~3.7 days is not interpretable. A tolerance below this
/// also runs into the branch-length grid, whose spacing quantizes node times -- movements come
/// out as integer multiples of a per-dataset step -- so a tighter criterion would make runs
/// exhaust `--max-iter` chasing jitter that no longer changes the answer.
pub const NODE_TIME_TOLERANCE_YEARS: f64 = 1e-2;

/// Tracks convergence metrics across timetree optimization iterations.
///
/// Records log-likelihood components and change counts to monitor convergence:
/// - Node time movement (max_time_change) should approach zero
/// - Polytomies resolved (n_resolved) should stabilize
/// - Log-likelihoods should increase or stabilize
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceMetrics {
  /// Number of ancestral sequence changes in this iteration
  pub n_diff: usize,
  /// Number of polytomies resolved in this iteration
  pub n_resolved: usize,
  /// Largest absolute change in any node's time over this iteration, in years.
  /// `None` when no node was dated both before and after the iteration.
  pub max_time_change: Option<f64>,
  /// Root-mean-square change in node times over this iteration, in years.
  /// Read against `max_time_change`: a large ratio means a few nodes are moving, not the tree.
  pub rms_time_change: Option<f64>,
  /// Sequence log-likelihood (log-probability of sequences given the tree and substitution model)
  pub log_lh_seq: Option<LogLh>,
  /// Positional log-likelihood (log-probability of node positions on the time axis)
  pub log_lh_pos: Option<LogLh>,
  /// Coalescent log-likelihood (population genetic log-prior on node times)
  pub log_lh_coal: Option<LogLh>,
  /// Total log-likelihood (sum of available components; absent components are excluded)
  pub log_lh_total: Option<LogLh>,
}

impl ConvergenceMetrics {
  /// Node times are what the refinement loop moves, so they are what convergence is judged on.
  ///
  /// `n_diff` counts changes in the reconstructed ancestral sequences, which is a coarser signal:
  /// the maximum-likelihood state at a site rarely flips even when every node time shifts, so a
  /// round can move the tree and still report `n_diff == 0`. It is kept as the fallback for the
  /// case where no node is dated on both sides of the iteration and no movement can be measured.
  pub fn has_converged(&self) -> bool {
    let times_settled = self
      .max_time_change
      .map_or(self.n_diff == 0, |change| change < NODE_TIME_TOLERANCE_YEARS);
    times_settled && self.n_resolved == 0
  }
}
