use crate::commands::optimize::optimize_dense;
use crate::commands::optimize::optimize_eval::evaluate_site_contributions;
use crate::commands::optimize::optimize_unified::OptimizationMetrics;

/// Evaluate dense contribution for a given branch length (with derivatives).
pub fn evaluate_dense_contribution(
  contribution: &optimize_dense::PartitionContribution,
  branch_length: f64,
) -> OptimizationMetrics {
  evaluate_dense_contribution_impl(contribution, branch_length, true)
}

/// Evaluate dense contribution for a given branch length (optional derivatives).
pub fn evaluate_dense_contribution_impl(
  contribution: &optimize_dense::PartitionContribution,
  branch_length: f64,
  compute_derivatives: bool,
) -> OptimizationMetrics {
  let sites = contribution.coefficients.outer_iter().map(|row| (1.0, row));
  evaluate_site_contributions(sites, &contribution.gtr.eigvals, branch_length, compute_derivatives)
}
