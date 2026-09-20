use crate::make_internal_error;
use crate::optimize::likelihood::{OptimizationMetrics, evaluate_with_indels};
use crate::optimize::zero_boundary::grid_search_inner;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use num::clamp;

pub(super) const NEWTON_REL_TOL: f64 = 0.001;

const NEWTON_ABS_TOL: f64 = 1e-8;

const NEWTON_MAX_ITER: usize = 10;

pub(super) fn newton_tolerance_t(t: f64) -> f64 {
  f64::max(NEWTON_REL_TOL * t, NEWTON_ABS_TOL)
}

pub(super) fn newton_tolerance_sqrt(s: f64) -> f64 {
  f64::max(0.5 * NEWTON_REL_TOL * s, NEWTON_ABS_TOL)
}

pub(super) fn newton_tolerance_log() -> f64 {
  NEWTON_REL_TOL.ln_1p()
}

pub(super) fn newton_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  if metrics.second_derivative < 0.0 {
    let mut bl = branch_length;
    let mut new_bl = (bl - clamp(metrics.derivative / metrics.second_derivative, -1.0, bl)).max(min_branch_length);

    for _ in 0..NEWTON_MAX_ITER {
      if (new_bl - bl).abs() <= newton_tolerance_t(bl) {
        break;
      }
      let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, new_bl)?;
      if new_metrics.second_derivative < 0.0 {
        bl = new_bl;
        new_bl = (bl - clamp(new_metrics.derivative / new_metrics.second_derivative, -1.0, bl)).max(min_branch_length);
      } else {
        break;
      }
    }
    Ok(new_bl)
  } else {
    grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation)
  }
}

pub(super) fn newton_sqrt_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  let mut s = branch_length.sqrt();
  let min_s = min_branch_length.sqrt();

  let (ds, d2s) = chain_rule_sqrt(s, metrics.derivative, metrics.second_derivative);

  if d2s >= 0.0 {
    return grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation);
  }

  let mut new_s = (s - clamp(ds / d2s, sqrt_step_lower_bound(s), s)).max(min_s);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_s - s).abs() <= newton_tolerance_sqrt(s) {
      break;
    }
    let t = new_s * new_s;
    let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, t)?;
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

pub(super) fn sqrt_step_lower_bound(s: f64) -> f64 {
  s - (s * s + 1.0).sqrt()
}

pub(super) fn log_step_lower_bound(t: f64) -> f64 {
  -(1.0 / t).ln_1p()
}

pub(super) fn chain_rule_sqrt(s: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_ds = 2.0 * s * dl_dt;
  let d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
  (dl_ds, d2l_ds2)
}

pub(super) fn chain_rule_log(t: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_du = t * dl_dt;
  let d2l_du2 = t * t * d2l_dt2 + t * dl_dt;
  (dl_du, d2l_du2)
}

pub(super) fn newton_log_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<f64, Report> {
  if !branch_length.is_finite() || branch_length <= 0.0 {
    return make_internal_error!(
      "newton_log_inner: branch_length must be finite and strictly positive, got {branch_length}"
    );
  }

  let mut u = branch_length.ln();
  let u_min = min_branch_length.max(1e-12).ln();

  let (du, d2u) = chain_rule_log(branch_length, metrics.derivative, metrics.second_derivative);

  if d2u >= 0.0 {
    return grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation);
  }

  let step_lower = log_step_lower_bound(branch_length);
  let step_upper = u - u_min;
  let mut new_u = (u - clamp(du / d2u, step_lower, step_upper)).max(u_min);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_u - u).abs() <= newton_tolerance_log() {
      break;
    }
    let t = new_u.exp();
    let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, t)?;
    let (du_new, d2u_new) = chain_rule_log(t, new_metrics.derivative, new_metrics.second_derivative);
    if d2u_new < 0.0 {
      u = new_u;
      let step_lower = log_step_lower_bound(t);
      let step_upper = u - u_min;
      new_u = (u - clamp(du_new / d2u_new, step_lower, step_upper)).max(u_min);
    } else {
      break;
    }
  }

  Ok(new_u.exp())
}
