use crate::commands::optimize::optimize_unified::OptimizationMetrics;
use ndarray::{Array1, ArrayView1};

/// Unified evaluation of eigenvalue-space site contributions.
///
/// Computes log-likelihood, first derivative, and second derivative of the
/// branch-length likelihood from a sequence of (multiplicity, coefficient)
/// pairs. Dense and sparse representations differ only in iteration: dense
/// yields `(1.0, row)` for each alignment position, sparse yields
/// `(multiplicity, coefficients)` for each compressed site pattern.
///
/// The per-site log-likelihood for branch length $t$ is:
///
///   $\ell_i(t) = \ln\!\bigl(\sum_c k_{ic}\, e^{\lambda_c t}\bigr)$
///
/// where $k_{ic}$ are the eigenvalue-space coefficients and $\lambda_c$ the
/// eigenvalues. The first and second derivatives follow from the quotient rule.
pub fn evaluate_site_contributions<'a>(
  sites: impl Iterator<Item = (f64, ArrayView1<'a, f64>)>,
  eigvals: &Array1<f64>,
  branch_length: f64,
  compute_derivatives: bool,
) -> OptimizationMetrics {
  let mut log_lh = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;

  let exp_ev = (eigvals * branch_length).mapv(f64::exp);

  if compute_derivatives {
    let ev_exp_ev = eigvals * &exp_ev;
    let ev2_exp_ev = eigvals * &ev_exp_ev;

    for (multiplicity, coefficients) in sites {
      let site_lh = (&coefficients * &exp_ev).sum();
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += multiplicity * site_lh.ln();
      let d1 = (&coefficients * &ev_exp_ev).sum() / site_lh;
      derivative += multiplicity * d1;
      second_derivative += multiplicity * (&coefficients * &ev2_exp_ev).sum() / site_lh - multiplicity * d1.powi(2);
    }
  } else {
    for (multiplicity, coefficients) in sites {
      let site_lh = (&coefficients * &exp_ev).sum();
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += multiplicity * site_lh.ln();
    }
  }

  OptimizationMetrics::new(log_lh, derivative, second_derivative)
}
