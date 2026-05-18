use crate::optimize::dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::optimize::indel::poisson_indel_log_lh;
use crate::optimize::sparse_eval::{evaluate_sparse_contribution, evaluate_sparse_contribution_impl};
use crate::partition::optimization_contribution::OptimizationContribution;

/// Metrics computed during branch length optimization
#[derive(Clone, Debug, Default)]
pub struct OptimizationMetrics {
  /// Log likelihood value
  pub log_lh: f64,
  /// First derivative (gradient) of log likelihood with respect to branch length
  pub derivative: f64,
  /// Second derivative (hessian) of log likelihood with respect to branch length
  pub second_derivative: f64,
}

impl OptimizationMetrics {
  pub fn new(log_lh: f64, derivative: f64, second_derivative: f64) -> Self {
    Self {
      log_lh,
      derivative,
      second_derivative,
    }
  }

  /// Add another set of metrics to this one
  pub fn add(&mut self, other: &OptimizationMetrics) {
    self.log_lh += other.log_lh;
    self.derivative += other.derivative;
    self.second_derivative += other.second_derivative;
  }
}

#[allow(clippy::multiple_inherent_impl)]
impl OptimizationContribution {
  pub fn evaluate(&self, branch_length: f64) -> OptimizationMetrics {
    match self {
      OptimizationContribution::Dense(contribution) => evaluate_dense_contribution(contribution, branch_length),
      OptimizationContribution::Sparse(contribution) => evaluate_sparse_contribution(contribution, branch_length),
    }
  }
}

pub fn evaluate_mixed(contributions: &[OptimizationContribution], branch_length: f64) -> OptimizationMetrics {
  evaluate_mixed_impl(contributions, branch_length, true)
}

pub fn evaluate_mixed_log_lh_only(contributions: &[OptimizationContribution], branch_length: f64) -> f64 {
  evaluate_mixed_impl(contributions, branch_length, false).log_lh
}

fn evaluate_mixed_impl(
  contributions: &[OptimizationContribution],
  branch_length: f64,
  compute_derivatives: bool,
) -> OptimizationMetrics {
  let mut total_metrics = OptimizationMetrics::default();
  for contribution in contributions {
    let metrics = match contribution {
      OptimizationContribution::Dense(c) => evaluate_dense_contribution_impl(c, branch_length, compute_derivatives),
      OptimizationContribution::Sparse(c) => evaluate_sparse_contribution_impl(c, branch_length, compute_derivatives),
    };
    total_metrics.add(&metrics);
  }
  total_metrics
}

/// Evaluate substitution + indel contributions for a given branch length.
pub fn evaluate_with_indels(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> OptimizationMetrics {
  let mut metrics = evaluate_mixed(contributions, branch_length);
  metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, branch_length));
  metrics
}

/// Log-likelihood only (no derivatives) for substitution + indel contributions.
pub fn evaluate_with_indels_log_lh_only(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> f64 {
  let sub_lh = evaluate_mixed_log_lh_only(contributions, branch_length);
  let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, branch_length).log_lh;
  sub_lh + indel_lh
}
