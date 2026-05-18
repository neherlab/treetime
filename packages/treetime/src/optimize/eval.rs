use crate::optimize::likelihood::OptimizationMetrics;
use itertools::izip;
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
/// eigenvalues. Define $S_i = \sum_c k_{ic} e^{\lambda_c t}$ (the site
/// likelihood) and the posterior weights $w_{ic} = k_{ic} e^{\lambda_c t} / S_i$.
/// The weights sum to 1 even when some $k_{ic}$ are negative, so
///
///   $\ell'_i(t)  = S'_i / S_i                 = \sum_c w_{ic}\, \lambda_c$
///   $\ell''_i(t) = S''_i / S_i - (S'_i/S_i)^2 = \sum_c w_{ic}\, (\lambda_c - \ell'_i)^2$
///
/// The Hessian is expressed as the posterior variance of the eigenvalues
/// rather than the difference $E[\lambda^2] - E[\lambda]^2$. The centered
/// form avoids the catastrophic cancellation that the difference suffers
/// when the posterior is tightly peaked (large $t$, degenerate spectra,
/// or sites close to stationarity), where both moments have nearly the
/// same magnitude.
#[allow(single_use_lifetimes)]
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

    for (multiplicity, coefficients) in sites {
      let k_exp = &coefficients * &exp_ev;
      let site_lh = k_exp.sum();
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += multiplicity * site_lh.ln();

      // First derivative: posterior mean eigenvalue, `S'/S`.
      let mean_ev = (&coefficients * &ev_exp_ev).sum() / site_lh;
      derivative += multiplicity * mean_ev;

      // Second derivative: posterior variance of eigenvalues in centered
      // (Welford) form, `sum_c w_c * (lambda_c - mean)^2`. Mathematically
      // equivalent to `S''/S - (S'/S)^2` but numerically stable: the
      // subtraction happens once between `lambda_c` and the mean (similar
      // absolute scale), not between two `O(lambda^2)` moments.
      let variance = izip!(k_exp.iter(), eigvals.iter())
        .map(|(&ke, &l)| {
          let d = l - mean_ev;
          ke * d * d
        })
        .sum::<f64>()
        / site_lh;
      second_derivative += multiplicity * variance;
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
