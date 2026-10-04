use crate::optimize::branch_length::validate_branch_length_value;
use crate::optimize::likelihood::OptimizationMetrics;
use eyre::Report;
use itertools::izip;
use ndarray::{Array1, ArrayView1, Zip};
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

  let mut k_exp = Array1::zeros(exp_ev.raw_dim());

  if compute_derivatives {
    let ev_exp_ev = eigvals * &exp_ev;
    let mut k_ev_exp = Array1::zeros(exp_ev.raw_dim());

    for (multiplicity, coefficients) in sites {
      multiply_into(&mut k_exp, &coefficients, &exp_ev);
      let site_lh = k_exp.sum();
      debug_assert!(site_lh.is_finite(), "Non-finite site likelihood: {site_lh}");
      log_lh += multiplicity * site_lh.ln();

      multiply_into(&mut k_ev_exp, &coefficients, &ev_exp_ev);
      let mean_ev = k_ev_exp.sum() / site_lh;
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
      multiply_into(&mut k_exp, &coefficients, &exp_ev);
      let site_lh = k_exp.sum();
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

fn multiply_into(out: &mut Array1<f64>, lhs: &ArrayView1<'_, f64>, rhs: &Array1<f64>) {
  Zip::from(out).and(lhs).and(rhs).for_each(|out, &lhs, &rhs| *out = lhs * rhs);
}
