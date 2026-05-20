use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::timetree::EdgeTimetree;
use crate::payload::timetree::NodeTimetree;
use crate::timetree::convergence::likelihood::{
  compute_coalescent_likelihood, compute_positional_likelihood, compute_sequence_likelihood,
};
use crate::timetree::convergence::metrics::ConvergenceMetrics;
use eyre::Report;
use log::info;
use parking_lot::RwLock;
use std::io::Write;
use std::sync::Arc;
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
    graph: &GraphTimetree,
    partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
    coalescent_tc: Option<&Distribution>,
  ) -> Result<(), Report> {
    let lh_seq = compute_sequence_likelihood(graph, partitions);
    let lh_pos = compute_positional_likelihood(graph);
    let lh_coal = compute_coalescent_likelihood(graph, coalescent_tc);
    let lh_total = [lh_seq, lh_pos, lh_coal].into_iter().flatten().reduce(|acc, v| acc + v);

    let metric = ConvergenceMetrics {
      n_diff,
      n_resolved,
      lh_seq,
      lh_pos,
      lh_coal,
      lh_total,
    };

    if let Some(writer) = &mut self.tracelog_writer {
      writer.write(&metric)?;
    }

    info!(
      "  Iteration {}: n_diff={n_diff}, n_resolved={n_resolved}, lh_seq={:.2}, lh_pos={:.2}, lh_coal={:.2}, total_LH={:.2}{}",
      self.i,
      metric.lh_seq.unwrap_or(f64::NAN),
      metric.lh_pos.unwrap_or(f64::NAN),
      metric.lh_coal.unwrap_or(f64::NAN),
      metric.lh_total.unwrap_or(f64::NAN),
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
