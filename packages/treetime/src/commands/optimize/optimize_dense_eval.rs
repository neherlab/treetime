use crate::commands::optimize::optimize_dense;
use crate::commands::optimize::optimize_unified::OptimizationMetrics;

/// Evaluate dense contribution for a given branch length
pub fn evaluate_dense_contribution(
  contribution: &optimize_dense::PartitionContribution,
  branch_length: f64,
) -> OptimizationMetrics {
  evaluate_dense_contribution_impl(contribution, branch_length, true)
}

/// Evaluate dense contribution for a given branch length (implementation with optional derivatives)
pub fn evaluate_dense_contribution_impl(
  contribution: &optimize_dense::PartitionContribution,
  branch_length: f64,
  compute_derivatives: bool,
) -> OptimizationMetrics {
  let mut log_lh = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;

  let gtr = &contribution.gtr;
  let coefficients = &contribution.coefficients;

  let exp_ev = gtr.exp_eigvals_branch_length(branch_length);

  if compute_derivatives {
    let ev_exp_ev = &gtr.eigvals * &exp_ev;
    let ev2_exp_ev = &gtr.eigvals * &ev_exp_ev;

    for coefficients in coefficients.outer_iter() {
      let site_lh = (&coefficients * &exp_ev).sum();
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += site_lh.ln();
      let d1 = (&coefficients * &ev_exp_ev).sum() / site_lh;
      derivative += d1;
      second_derivative += (&coefficients * &ev2_exp_ev).sum() / site_lh - d1.powi(2);
    }
  } else {
    for coefficients in coefficients.outer_iter() {
      let site_lh = (&coefficients * &exp_ev).sum();
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += site_lh.ln();
    }
  }

  OptimizationMetrics::new(log_lh, derivative, second_derivative)
}
