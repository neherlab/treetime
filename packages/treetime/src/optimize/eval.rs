use crate::optimize::branch_length::validate_branch_length_value;
use crate::optimize::likelihood::OptimizationMetrics;
use eyre::Report;
use itertools::izip;
use ndarray::{Array1, ArrayView1};
use treetime_primitives::LogLh;

pub(crate) fn evaluate_site_contributions<'a>(
  sites: impl Iterator<Item = (f64, ArrayView1<'a, f64>)>,
  rate_eigvals: &Array1<f64>,
  branch_length: f64,
  compute_derivatives: bool,
) -> Result<OptimizationMetrics, Report> {
  validate_branch_length_value(branch_length)?;

  let mut log_lh = 0.0;
  let mut derivative = 0.0;
  let mut second_derivative = 0.0;

  let exp_ev = rate_eigvals.mapv(|l| (l * branch_length).exp());

  if compute_derivatives {
    let ev_exp_ev = rate_eigvals * &exp_ev;

    for (multiplicity, coefficients) in sites {
      let site_lh = site_likelihood(&coefficients, &exp_ev);
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += multiplicity * site_lh.ln();

      let mean_ev = coefficients.dot(&ev_exp_ev) / site_lh;
      derivative += multiplicity * mean_ev;

      let variance = izip!(&coefficients, &exp_ev, rate_eigvals)
        .map(|(&c, &e, &l)| {
          let d = l - mean_ev;
          c * e * d * d
        })
        .sum::<f64>()
        / site_lh;
      second_derivative += multiplicity * variance;
    }
  } else {
    for (multiplicity, coefficients) in sites {
      let site_lh = site_likelihood(&coefficients, &exp_ev);
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += multiplicity * site_lh.ln();
    }
  }

  Ok(OptimizationMetrics::new(
    LogLh::new(log_lh),
    derivative,
    second_derivative,
  ))
}

fn site_likelihood(coefficients: &ArrayView1<'_, f64>, exp_ev: &Array1<f64>) -> f64 {
  coefficients.dot(exp_ev).max(f64::MIN_POSITIVE)
}
