use crate::branch_lengths::branch_lengths_or_zero;
use crate::optimize::indel::EdgeIndels;
use crate::partition::set::PartitionSet;
use eyre::Report;
use log::debug;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_primitives::LogLh;
use treetime_utils::fmt::float::float_to_significant_digits;

pub(crate) struct TopologyHistory {
  dp: f64,
  log_lhs: Vec<LogLh>,
  best: Option<BestIteration>,
}

impl TopologyHistory {
  pub(crate) const fn new(dp: f64) -> Self {
    Self {
      dp,
      log_lhs: Vec::new(),
      best: None,
    }
  }

  pub(crate) fn record(
    &mut self,
    iteration: usize,
    log_lh: LogLh,
    branch_lengths: &EdgeMap<Option<f64>>,
    partitions: &PartitionSet,
    indels: &EdgeIndels,
  ) -> Option<ConvergenceReason> {
    self.log_lhs.push(log_lh);
    if !log_lh.value().is_finite() {
      return Some(ConvergenceReason::NumericalFailure);
    }
    if self.best.as_ref().is_none_or(|best| log_lh > best.log_lh) {
      self.best = Some(BestIteration {
        iteration,
        log_lh,
        branch_lengths: branch_lengths.clone(),
        rates: partitions.rates(),
        indels: indels.clone(),
      });
    }
    let best_lh = self.best.as_ref().map_or(log_lh, |best| best.log_lh);
    stop_reason(&self.log_lhs, best_lh, self.dp)
  }

  pub(crate) fn topology_changed(&mut self) {
    self.best = None;
  }

  pub(crate) fn into_best(self) -> Option<BestIteration> {
    self.best
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ConvergenceReason {
  Converged,
  Oscillating,
  Worsened,
  NumericalFailure,
}

pub(crate) struct BestIteration {
  iteration: usize,
  log_lh: LogLh,
  branch_lengths: EdgeMap<Option<f64>>,
  rates: Vec<f64>,
  indels: EdgeIndels,
}

impl BestIteration {
  pub(crate) const fn iteration(&self) -> usize {
    self.iteration
  }

  pub(crate) fn restore(
    self,
    graph: &Graph,
    mut partitions: PartitionSet,
    branch_lengths: &mut EdgeMap<Option<f64>>,
    indels: &mut EdgeIndels,
  ) -> Result<PartitionSet, Report> {
    debug!(
      "Returning the branch lengths and rates of iteration {} with the highest likelihood {}",
      self.iteration + 1,
      float_to_significant_digits(self.log_lh.value(), 7)
    );
    *branch_lengths = self.branch_lengths;
    partitions.set_rates(&self.rates)?;
    *indels = self.indels;
    Ok(
      partitions
        .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
        .0,
    )
  }
}

fn stop_reason(lh_history: &[LogLh], best_lh: LogLh, dp: f64) -> Option<ConvergenceReason> {
  let iteration = lh_history.len().checked_sub(1)?;
  let earlier = |steps: usize| {
    iteration
      .checked_sub(steps)
      .map_or(LogLh::IMPOSSIBLE, |i| lh_history[i])
  };
  let lh = lh_history[iteration];
  let lh_prev = earlier(1);
  let lh_prev_prev = earlier(2);
  if (lh - lh_prev).abs() < dp.abs() {
    Some(ConvergenceReason::Converged)
  } else if iteration >= 2 && (lh - lh_prev_prev).abs() < dp.abs() {
    Some(ConvergenceReason::Oscillating)
  } else if iteration >= 2 && lh < lh_prev && lh_prev >= best_lh {
    Some(ConvergenceReason::Worsened)
  } else {
    None
  }
}
