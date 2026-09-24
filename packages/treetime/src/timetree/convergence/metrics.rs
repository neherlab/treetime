use serde::{Deserialize, Serialize};
use treetime_primitives::LogLh;

pub(crate) const NODE_TIME_TOLERANCE_YEARS: f64 = 1e-2;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceMetrics {
  pub(crate) n_diff: usize,
  pub(crate) n_resolved: usize,
  pub(crate) max_time_change: Option<f64>,
  pub(crate) rms_time_change: Option<f64>,
  pub(crate) log_lh_seq: Option<LogLh>,
  pub(crate) log_lh_pos: Option<LogLh>,
  pub(crate) log_lh_coal: Option<LogLh>,
  pub(crate) log_lh_total: Option<LogLh>,
}

impl ConvergenceMetrics {
  pub(crate) fn has_converged(&self) -> bool {
    let times_settled = self
      .max_time_change
      .map_or(self.n_diff == 0, |change| change < NODE_TIME_TOLERANCE_YEARS);
    times_settled && self.n_resolved == 0
  }
}
