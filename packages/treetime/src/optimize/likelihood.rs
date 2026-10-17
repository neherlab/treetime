use crate::optimize::dense_eval::evaluate_dense_contribution;
use crate::optimize::indel::IndelTerm;
use crate::optimize::sparse_eval::evaluate_sparse_contribution;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use treetime_primitives::LogLh;

pub(crate) fn evaluate_with_indels(
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  let mut metrics = evaluate_partitions(contributions, branch_length)?;
  metrics.add(&indels.metrics(branch_length)?);
  Ok(metrics)
}

pub(crate) fn evaluate_partitions(
  contributions: &[OptimizationContribution],
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  sum_partition_metrics(contributions, branch_length, true)
}

pub(crate) fn evaluate_with_indels_log_lh_only(
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  branch_length: f64,
) -> Result<LogLh, Report> {
  let sub_lh = evaluate_partitions_log_lh_only(contributions, branch_length)?;
  let indel_lh = indels.metrics(branch_length)?.log_lh;
  Ok(sub_lh + indel_lh)
}

pub(crate) fn evaluate_partitions_log_lh_only(
  contributions: &[OptimizationContribution],
  branch_length: f64,
) -> Result<LogLh, Report> {
  Ok(sum_partition_metrics(contributions, branch_length, false)?.log_lh)
}

fn sum_partition_metrics(
  contributions: &[OptimizationContribution],
  branch_length: f64,
  compute_derivatives: bool,
) -> Result<OptimizationMetrics, Report> {
  let mut total_metrics = OptimizationMetrics::default();
  for contribution in contributions {
    let metrics = match contribution {
      OptimizationContribution::Dense(c) => evaluate_dense_contribution(c, branch_length, compute_derivatives),
      OptimizationContribution::Sparse(c) => evaluate_sparse_contribution(c, branch_length, compute_derivatives),
    }?;
    total_metrics.add(&metrics);
  }
  Ok(total_metrics)
}

#[derive(Clone, Debug, Default)]
pub struct OptimizationMetrics {
  pub(crate) log_lh: LogLh,
  pub(crate) derivative: f64,
  pub(crate) second_derivative: f64,
}

impl OptimizationMetrics {
  pub(crate) fn new(log_lh: LogLh, derivative: f64, second_derivative: f64) -> Self {
    Self {
      log_lh,
      derivative,
      second_derivative,
    }
  }

  pub(crate) fn add(&mut self, other: &OptimizationMetrics) {
    self.log_lh += other.log_lh;
    self.derivative += other.derivative;
    self.second_derivative += other.second_derivative;
  }
}
