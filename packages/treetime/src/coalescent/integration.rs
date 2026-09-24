use eyre::Report;
use ndarray::Array1;
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;
use treetime_utils::make_error;

pub(super) fn compute_merger_rate_total_scalar(k: f64, tc: f64) -> f64 {
  compute_merger_rate_total(compute_merger_rate_lineage_count(k), tc)
}

pub(super) fn compute_integral_merger_rate(
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
) -> Result<PiecewiseLinearFn, Report> {
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.len() < 2 {
    return make_error!("lineage count must have at least 2 breakpoints");
  }

  let n = breakpoints.len();
  let mut integral_values = vec![0.0; n];

  for i in (0..n - 1).rev() {
    let t0 = breakpoints[i];
    let t1 = breakpoints[i + 1];
    let dt = t1 - t0;

    let mid = f64::midpoint(t0, t1);
    let k = lineage_counts.eval(mid);
    let tc = tc_dist.eval(mid)?;

    if !k.is_finite() {
      return make_error!("Coalescent lineage count must be finite at calendar time {mid:.6e}, got {k:.6e}");
    }
    if !tc.is_finite() || tc <= 0.0 {
      return make_error!("Coalescent Tc must be finite and positive at calendar time {mid:.6e}, got {tc:.6e}");
    }
    let rate = compute_merger_rate_per_lineage_scalar(k, tc);
    integral_values[i] = integral_values[i + 1] + dt * rate;
  }

  Ok(PiecewiseLinearFn::new(
    Array1::from(breakpoints.to_vec()),
    Array1::from_vec(integral_values),
  ))
}

pub(super) fn compute_merger_rate_per_lineage_scalar(k: f64, tc: f64) -> f64 {
  compute_merger_rate_per_lineage(compute_merger_rate_lineage_count(k), tc)
}

fn compute_merger_rate_lineage_count(k: f64) -> f64 {
  let n_merge_candidates = k - 1.0;
  if n_merge_candidates.is_nan() {
    n_merge_candidates
  } else {
    f64::max(0.5, n_merge_candidates)
  }
}

fn compute_merger_rate_per_lineage(n_merge_candidates: f64, tc: f64) -> f64 {
  0.5 * n_merge_candidates / tc
}

fn compute_merger_rate_total(n_merge_candidates: f64, tc: f64) -> f64 {
  0.5 * n_merge_candidates * (n_merge_candidates + 1.0) / tc
}
