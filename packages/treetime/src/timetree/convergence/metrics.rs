use crate::clock::clock_model::ClockModel;
use deser::{Deserialize, Serialize};
use treetime_primitives::LogLh;

pub(crate) const NODE_TIME_TOLERANCE_YEARS: f64 = 1e-2;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceMetrics {
  pub n_diff: usize,
  pub n_resolved: usize,
  pub max_time_change: Option<f64>,
  pub rms_time_change: Option<f64>,
  pub log_lh_seq: Option<LogLh>,
  pub log_lh_pos: Option<LogLh>,
  pub log_lh_coal: Option<LogLh>,
  pub log_lh_total: Option<LogLh>,
}

impl ConvergenceMetrics {
  pub(crate) fn has_converged(&self) -> bool {
    let times_settled = self
      .max_time_change
      .map_or(self.n_diff == 0, |change| change < NODE_TIME_TOLERANCE_YEARS);
    times_settled && self.n_resolved == 0
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct IterationClock {
  pub clock_rate: f64,
  pub r_squared: Option<f64>,
}

impl IterationClock {
  pub(crate) fn of(clock_model: &ClockModel) -> Self {
    Self {
      clock_rate: clock_model.clock_rate(),
      r_squared: clock_model.r_squared(),
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IterationRecord {
  pub iteration: usize,
  pub metrics: ConvergenceMetrics,
  pub clock: IterationClock,
}
