use crate::optimize::eval::evaluate_site_contributions;
use crate::optimize::likelihood::OptimizationMetrics;
use crate::partition::optimize_sparse;
use eyre::Report;

/// Evaluate sparse contribution for a given branch length (with derivatives).
pub fn evaluate_sparse_contribution(
  contribution: &optimize_sparse::PartitionContribution,
  branch_length: f64,
) -> Result<OptimizationMetrics, Report> {
  evaluate_sparse_contribution_impl(contribution, branch_length, true)
}

/// Evaluate sparse contribution for a given branch length (optional derivatives).
pub fn evaluate_sparse_contribution_impl(
  contribution: &optimize_sparse::PartitionContribution,
  branch_length: f64,
  compute_derivatives: bool,
) -> Result<OptimizationMetrics, Report> {
  let sites = contribution
    .site_contributions
    .iter()
    .map(|sc| (sc.multiplicity, sc.coefficients.view()));
  evaluate_site_contributions(sites, &contribution.gtr.eigvals, branch_length, compute_derivatives)
}
