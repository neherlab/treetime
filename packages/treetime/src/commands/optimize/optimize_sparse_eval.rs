use crate::commands::optimize::optimize_sparse;
use crate::commands::optimize::optimize_unified::OptimizationMetrics;

/// Evaluate sparse contribution for a given branch length
pub fn evaluate_sparse_contribution(
  contribution: &optimize_sparse::PartitionContribution,
  branch_length: f64,
) -> OptimizationMetrics {
  evaluate_sparse_contribution_impl(contribution, branch_length, true)
}

/// Evaluate sparse contribution for a given branch length (implementation with optional derivatives)
pub fn evaluate_sparse_contribution_impl(
  contribution: &optimize_sparse::PartitionContribution,
  branch_length: f64,
  compute_derivatives: bool,
) -> OptimizationMetrics {
  let mut log_lh = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;

  let exp_ev = contribution.gtr.exp_eigvals_branch_length(branch_length);

  if compute_derivatives {
    let ev_exp_ev = &contribution.gtr.eigvals * &exp_ev;
    let ev2_exp_ev = &contribution.gtr.eigvals * &ev_exp_ev;

    for optimize_sparse::SiteContribution {
      multiplicity,
      coefficients,
    } in &contribution.site_contributions
    {
      let site_lh = (coefficients * &exp_ev).sum();
      log_lh += multiplicity * site_lh.ln();
      let d1 = (coefficients * &ev_exp_ev).sum() / site_lh;
      derivative += multiplicity * d1;
      second_derivative += multiplicity * (coefficients * &ev2_exp_ev).sum() / site_lh - multiplicity * d1.powi(2);
    }
  } else {
    for optimize_sparse::SiteContribution {
      multiplicity,
      coefficients,
    } in &contribution.site_contributions
    {
      let site_lh = (coefficients * &exp_ev).sum();
      log_lh += multiplicity * site_lh.ln();
    }
  }

  OptimizationMetrics::new(log_lh, derivative, second_derivative)
}
