use crate::optimize::dense_eval::{evaluate_dense_contribution, evaluate_dense_contribution_impl};
use crate::optimize::indel::poisson_indel_log_lh;
use crate::optimize::sparse_eval::{evaluate_sparse_contribution, evaluate_sparse_contribution_impl};
use crate::partition::optimization_contribution::OptimizationContribution;
use eyre::Report;
use treetime_primitives::LogLh;

/// Metrics computed during branch length optimization
#[derive(Clone, Debug, Default)]
pub struct OptimizationMetrics {
  /// Log likelihood value
  pub log_lh: LogLh,
  /// First derivative (gradient) of log likelihood with respect to branch length
  pub derivative: f64,
  /// Second derivative (hessian) of log likelihood with respect to branch length
  pub second_derivative: f64,
}

impl OptimizationMetrics {
  pub fn new(log_lh: LogLh, derivative: f64, second_derivative: f64) -> Self {
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
  pub fn evaluate(&self, branch_length: f64) -> Result<OptimizationMetrics, Report> {
    match self {
      OptimizationContribution::Dense(contribution) => evaluate_dense_contribution(contribution, branch_length),
      OptimizationContribution::Sparse(contribution) => evaluate_sparse_contribution(contribution, branch_length),
    }
  }
}

pub fn evaluate_mixed(
  contributions: &[OptimizationContribution],
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  evaluate_mixed_impl(contributions, branch_length, true)
}

pub fn evaluate_mixed_log_lh_only(
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

/// Evaluate substitution + indel contributions for a given branch length.
pub fn evaluate_with_indels(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  let mut metrics = evaluate_mixed(contributions, branch_length)?;
  metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, branch_length)?);
  Ok(metrics)
}

/// Log-likelihood only (no derivatives) for substitution + indel contributions.
pub fn evaluate_with_indels_log_lh_only(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  branch_length: f64,
) -> Result<LogLh, Report> {
  let sub_lh = evaluate_mixed_log_lh_only(contributions, branch_length)?;
  let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, branch_length)?.log_lh;
  Ok(sub_lh + indel_lh)
}

/// Boundary neg-log ordinate of the branch-length likelihood at `t = 0`, peak-normalized against
/// `max_log_lh`.
///
/// Classifies the hard boundary from the likelihood model: `Some(max_log_lh - log_lh(0))` when the
/// density is finite there (a zero-mutation branch, whose neg-log approaches the boundary as a
/// straight line), and `None` when it diverges (a forced substitution or an indel, whose neg-log
/// approaches as a power law). The returned value is the input a straight-line
/// [`HardApproachLaw`](treetime_grid::HardApproachLaw) fit needs; the `None` case selects the
/// power-law fit instead.
///
/// An indel branch cannot be evaluated at `t = 0` (`k*ln(0) = -inf`, and `poisson_indel_log_lh`
/// rejects `t = 0` for `k > 0`), so `indel_count > 0` diverges without evaluation.
pub fn branch_length_boundary_ordinate(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  max_log_lh: f64,
) -> Result<Option<f64>, Report> {
  if indel_count > 0 {
    return Ok(None);
  }
  let boundary_log_lh = evaluate_with_indels_log_lh_only(contributions, 0, indel_rate, 0.0)?.value();
  Ok(boundary_log_lh.is_finite().then_some(max_log_lh - boundary_log_lh))
}
