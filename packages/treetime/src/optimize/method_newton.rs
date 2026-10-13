use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels;
use crate::optimize::zero_boundary::grid_search_inner;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use num::clamp;

const NEWTON_REL_TOL: f64 = 0.001;

const NEWTON_ABS_TOL: f64 = 1e-8;

const NEWTON_MAX_ITER: usize = 10;

pub(super) fn newton_sqrt_inner(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let metrics = evaluate_with_indels(contributions, indels, branch_length)?;
  let mut s = branch_length.sqrt();
  let min_s = min_branch_length.sqrt();

  let (ds, d2s) = chain_rule_sqrt(s, metrics.derivative, metrics.second_derivative);

  if d2s >= 0.0 {
    return grid_search_inner(branch_length, contributions, indels, one_mutation);
  }

  let mut new_s = (s - clamp(ds / d2s, sqrt_step_lower_bound(s), s)).max(min_s);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_s - s).abs() <= newton_tolerance_sqrt(s) {
      break;
    }
    let t = new_s * new_s;
    let new_metrics = evaluate_with_indels(contributions, indels, t)?;
    let (ds_new, d2s_new) = chain_rule_sqrt(new_s, new_metrics.derivative, new_metrics.second_derivative);
    if d2s_new < 0.0 {
      s = new_s;
      new_s = (s - clamp(ds_new / d2s_new, sqrt_step_lower_bound(s), s)).max(min_s);
    } else {
      break;
    }
  }

  Ok(new_s * new_s)
}

pub(super) fn newton_tolerance_sqrt(s: f64) -> f64 {
  f64::max(0.5 * NEWTON_REL_TOL * s, NEWTON_ABS_TOL)
}

pub(super) fn sqrt_step_lower_bound(s: f64) -> f64 {
  s - (s * s + 1.0).sqrt()
}

pub(super) fn chain_rule_sqrt(s: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_ds = 2.0 * s * dl_dt;
  let d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
  (dl_ds, d2l_ds2)
}
