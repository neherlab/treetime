use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetreeRef};
use crate::timetree::convergence::likelihood::{
  compute_coalescent_log_lh, compute_positional_log_lh, compute_sequence_log_lh,
};
use crate::timetree::convergence::metrics::ConvergenceMetrics;
use crate::timetree::convergence::node_times::NodeTimeChange;
use eyre::Report;
use log::info;
use std::io::Write;
use treetime_distribution::Distribution;
use treetime_io::csv::CsvStructWriter;

pub struct TimetreeOptimizer {
  trace: Vec<ConvergenceMetrics>,
  tracelog_writer: Option<CsvStructWriter<Box<dyn Write + Send>>>,
  max_iterations: usize,
  suppress_convergence: bool,
  i: usize,
}

impl TimetreeOptimizer {
  pub fn new(max_iter: usize, suppress_convergence: bool) -> Self {
    Self {
      trace: vec![],
      tracelog_writer: None,
      max_iterations: max_iter,
      suppress_convergence,
      i: 0,
    }
  }

  pub fn with_tracelog(mut self, writer: impl Write + Send + 'static) -> Result<Self, Report> {
    let boxed: Box<dyn Write + Send> = Box::new(writer);
    self.tracelog_writer = Some(CsvStructWriter::new(boxed, b',')?);
    Ok(self)
  }

  pub fn next_iter(&mut self) -> Option<IterationContext> {
    if self.has_converged() || self.has_reached_max_iterations() {
      return None;
    }

    self.i += 1;
    info!("### Timetree iteration {}/{}", self.i, self.max_iterations);

    Some(IterationContext { i: self.i })
  }

  pub fn record(
    &mut self,
    n_diff: usize,
    n_resolved: usize,
    time_change: NodeTimeChange,
    graph: &GraphTimetree,
    partitions: &[PartitionTimetreeRef],
    coalescent_tc: Option<&Distribution>,
  ) -> Result<(), Report> {
    let log_lh_seq = compute_sequence_log_lh(graph, partitions);
    let log_lh_pos = compute_positional_log_lh(graph);
    let log_lh_coal = compute_coalescent_log_lh(graph, coalescent_tc);
    let log_lh_total = [log_lh_seq, log_lh_pos, log_lh_coal]
      .into_iter()
      .flatten()
      .reduce(|acc, v| acc + v);

    let metric = ConvergenceMetrics {
      n_diff,
      n_resolved,
      max_time_change: time_change.max,
      rms_time_change: time_change.rms,
      log_lh_seq,
      log_lh_pos,
      log_lh_coal,
      log_lh_total,
    };

    if let Some(writer) = &mut self.tracelog_writer {
      writer.write(&metric)?;
    }

    info!(
      "  Iteration {}: max_dt={:.4}, rms_dt={:.4}, n_diff={n_diff}, n_resolved={n_resolved}, log_lh_seq={:.2}, log_lh_pos={:.2}, log_lh_coal={:.2}, log_lh_total={:.2}{}",
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

  pub fn iteration_count(&self) -> usize {
    self.i
  }

  pub fn trace(&self) -> &[ConvergenceMetrics] {
    &self.trace
  }

  fn has_converged(&self) -> bool {
    !self.suppress_convergence && self.trace.last().is_some_and(|m| m.has_converged())
  }

  fn has_reached_max_iterations(&self) -> bool {
    self.i >= self.max_iterations
  }
}

pub struct IterationContext {
  pub i: usize,
}
