use crate::optimize::eval::evaluate_site_contributions;
use crate::optimize::likelihood::OptimizationMetrics;
use crate::partition::optimize;
use eyre::Report;

pub(crate) fn evaluate_sparse_contribution(
  contribution: &optimize::sparse::PartitionContribution,
  branch_length: f64,
  compute_derivatives: bool,
) -> Result<OptimizationMetrics, Report> {
  let sites = contribution
    .site_contributions
    .iter()
    .map(|sc| (sc.multiplicity, sc.coefficients.view()));
  evaluate_site_contributions(sites, &contribution.gtr.eigvals, branch_length, compute_derivatives)
}
