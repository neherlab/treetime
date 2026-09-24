use crate::optimize::branch_length::validate_branch_length_value;
use crate::optimize::likelihood::OptimizationMetrics;
use eyre::Report;
use itertools::izip;
use ndarray::{Array1, ArrayView1};
use treetime_primitives::LogLh;

#[allow(
  single_use_lifetimes,
  reason = "stable Rust cannot elide a lifetime inside an argument-position impl Trait"
)]
pub(crate) fn evaluate_site_contributions<'a>(
  sites: impl Iterator<Item = (f64, ArrayView1<'a, f64>)>,
  eigvals: &Array1<f64>,
  branch_length: f64,
  compute_derivatives: bool,
) -> Result<OptimizationMetrics, Report> {
  validate_branch_length_value(branch_length)?;

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

      let mean_ev = (&coefficients * &ev_exp_ev).sum() / site_lh;
      derivative += multiplicity * mean_ev;

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

  Ok(OptimizationMetrics::new(
    LogLh::new(log_lh),
    derivative,
    second_derivative,
  ))
}
