use crate::optimize::dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::optimize::indel::poisson_indel_log_lh;
use crate::optimize::sparse_eval::{evaluate_sparse_contribution, evaluate_sparse_contribution_impl};
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use treetime_primitives::LogLh;

#[derive(Clone, Debug, Default)]
pub struct OptimizationMetrics {
  pub log_lh: LogLh,
  pub derivative: f64,
  pub second_derivative: f64,
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

#[allow(clippy::multiple_inherent_impl)]
impl OptimizationContribution {
  pub fn evaluate(&self, branch_length: f64) -> Result<OptimizationMetrics, Report> {
    match self {
      OptimizationContribution::Dense(contribution) => evaluate_dense_contribution(contribution, branch_length),
      OptimizationContribution::Sparse(contribution) => evaluate_sparse_contribution(contribution, branch_length),
    }
  }
}

pub(crate) fn evaluate_mixed(
  contributions: &[OptimizationContribution],
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  evaluate_mixed_impl(contributions, branch_length, true)
}

pub(crate) fn evaluate_mixed_log_lh_only(
  contributions: &[OptimizationContribution],
  branch_length: f64,
) -> Result<LogLh, Report> {
  Ok(evaluate_mixed_impl(contributions, branch_length, false)?.log_lh)
}

fn evaluate_mixed_impl(
  contributions: &[OptimizationContribution],
  branch_length: f64,
  compute_derivatives: bool,
) -> Result<OptimizationMetrics, Report> {
  let mut total_metrics = OptimizationMetrics::default();
  for contribution in contributions {
    let metrics = match contribution {
      OptimizationContribution::Dense(c) => evaluate_dense_contribution_impl(c, branch_length, compute_derivatives),
      OptimizationContribution::Sparse(c) => evaluate_sparse_contribution_impl(c, branch_length, compute_derivatives),
    }?;
    total_metrics.add(&metrics);
  }
  Ok(total_metrics)
}

pub(crate) fn evaluate_with_indels(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  let mut metrics = evaluate_mixed(contributions, branch_length)?;
  metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, branch_length)?);
  Ok(metrics)
}

pub(crate) fn evaluate_with_indels_log_lh_only(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> Result<LogLh, Report> {
  let sub_lh = evaluate_mixed_log_lh_only(contributions, branch_length)?;
  let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, branch_length)?.log_lh;
  Ok(sub_lh + indel_lh)
}
