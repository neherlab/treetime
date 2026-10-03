use crate::progress::LogSink;
use crate::progress_info;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::convergence::likelihood::{
  compute_coalescent_log_lh, compute_positional_log_lh, compute_sequence_log_lh,
};
use crate::timetree::convergence::metrics::{ConvergenceMetrics, IterationClock, IterationRecord};
use crate::timetree::convergence::node_times::NodeTimeChange;
use crate::timetree::inference::result::TimeInference;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) struct TimetreeOptimizer<'a> {
  pub(crate) trace: Vec<ConvergenceMetrics>,
  trace_sink: Option<&'a mut dyn TraceSink>,
  max_iterations: usize,
  pub(crate) i: usize,
}

impl<'a> TimetreeOptimizer<'a> {
  pub(crate) fn new(max_iter: usize) -> Self {
    Self {
      trace: vec![],
      trace_sink: None,
      max_iterations: max_iter,
      i: 0,
    }
  }

  #[must_use]
  pub(crate) fn with_trace_sink(mut self, sink: &'a mut dyn TraceSink) -> Self {
    self.trace_sink = Some(sink);
    self
  }

  pub(crate) fn next_iter(&mut self, log: &dyn LogSink) -> Option<IterationContext> {
    if self.has_converged() || self.has_reached_max_iterations() {
      return None;
    }

    self.i += 1;
    progress_info!(log, "### Timetree iteration {}/{}", self.i, self.max_iterations);

    Some(IterationContext { i: self.i })
  }

  pub(crate) fn record(
    &mut self,
    sequence_changes: usize,
    resolved_nodes: usize,
    time_change: NodeTimeChange,
    graph: &Graph,
    branch_model: &BranchModel,
    inference: &TimeInference,
    coalescent_tc: Option<&Distribution>,
    clock: IterationClock,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    log: &dyn LogSink,
  ) -> Result<(), Report> {
    let log_lh_seq = compute_sequence_log_lh(graph, branch_model);
    let log_lh_pos = compute_positional_log_lh(graph, inference);
    let log_lh_coal =
      coalescent_tc.and_then(|tc| compute_coalescent_log_lh(graph, tc, &inference.coalescent_node_times(), names, log));
    let log_lh_total = [log_lh_seq, log_lh_pos, log_lh_coal]
      .into_iter()
      .flatten()
      .reduce(|acc, v| acc + v);

    let metric = ConvergenceMetrics {
      n_diff: sequence_changes,
      n_resolved: resolved_nodes,
      max_time_change: time_change.max,
      rms_time_change: time_change.rms,
      log_lh_seq,
      log_lh_pos,
      log_lh_coal,
      log_lh_total,
    };

    if let Some(sink) = &mut self.trace_sink {
      sink.emit(&IterationRecord {
        iteration: self.i,
        metrics: metric.clone(),
        clock,
      })?;
    }

    progress_info!(
      log,
      "  Iteration {}: max_dt={:.4}, rms_dt={:.4}, n_diff={sequence_changes}, n_resolved={resolved_nodes}, log_lh_seq={:.2}, log_lh_pos={:.2}, log_lh_coal={:.2}, log_lh_total={:.2}{}",
      self.i,
      metric.max_time_change.unwrap_or(f64::NAN),
      metric.rms_time_change.unwrap_or(f64::NAN),
      metric.log_lh_seq.map_or(f64::NAN, |log_lh| log_lh.value()),
      metric.log_lh_pos.map_or(f64::NAN, |log_lh| log_lh.value()),
      metric.log_lh_coal.map_or(f64::NAN, |log_lh| log_lh.value()),
      metric.log_lh_total.map_or(f64::NAN, |log_lh| log_lh.value()),
      if metric.has_converged() { " [converged]" } else { "" }
    );

    self.trace.push(metric);
    Ok(())
  }

  fn has_converged(&self) -> bool {
    self.trace.last().is_some_and(|m| m.has_converged())
  }

  fn has_reached_max_iterations(&self) -> bool {
    self.i >= self.max_iterations
  }
}

pub trait TraceSink: Send {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report>;
}

pub(crate) struct IterationContext {
  pub i: usize,
}
