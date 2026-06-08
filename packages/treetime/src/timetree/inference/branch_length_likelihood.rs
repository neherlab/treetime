use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimization_contribution::OptimizationContribution;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;

/// Compute the branch-length likelihood distribution used for time inference.
///
/// The grid likelihood includes both the substitution contribution (from
/// `contributions`) and the Poisson indel contribution (from `indel_count` and
/// `indel_rate`). For datasets without indels (`indel_rate == 0`), the Poisson
/// term is a no-op and the distribution reduces to the substitution likelihood.
///
/// Keeping the timetree grid consistent with the per-edge Newton evaluator in
/// `run_optimize_mixed()` is required for branches whose only evolutionary
/// signal is an indel event: otherwise the distribution peaks at zero while
/// the optimizer correctly assigns positive length (see
/// `docs/algorithms/optimize.md` and
/// `docs/port-intentional-changes/optimize-indel-contribution-to-likelihood.md`).
pub fn compute_branch_length_distribution(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  current_branch_length: f64,
  one_mutation: f64,
  n_grid_points: usize,
  clock_rate: f64,
  gamma: f64,
) -> Result<Arc<Distribution>, Report> {
  debug_assert!(clock_rate > 0.0, "clock_rate must be positive, got {clock_rate:.6e}");
  debug_assert!(gamma > 0.0);

  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points, clock_rate);

  // `create_simple_grid` always returns strictly positive branch lengths
  // (`min_bl = one_mutation * 0.1 > 0`), satisfying the `t > 0` precondition
  // of `poisson_indel_log_lh` when `indel_count > 0`.
  let log_lh =
    grid.mapv(|branch_len| evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, branch_len));
  let max_log_lh = log_lh.max()?;

  let normalized_prob = (&log_lh - *max_log_lh).exp();

  // Convert branch length grid to time grid: time = branch_length / (clock_rate * gamma)
  // gamma > 1 means faster evolution, so same substitutions correspond to shorter time
  let effective_clock_rate = clock_rate * gamma;
  let time_min = grid[0] / effective_clock_rate;
  let time_max = grid[grid.len() - 1] / effective_clock_rate;

  let distribution_fn = DistributionFunction::from_range_values((time_min, time_max), normalized_prob)?;
  Ok(Arc::new(Distribution::Function(distribution_fn)))
}

/// Maximum branch time in years for the distribution grid.
///
/// Controls how far the branch length grid extends by converting to subs/site
/// via clock_rate: max_bl = MAX_BRANCH_TIME * clock_rate. This is clock-rate
/// adaptive: fast-evolving datasets (flu, rate ~0.003) get max_bl ≈ 0.6;
/// slow-evolving datasets (ebola, rate ~0.0006) get max_bl ≈ 0.125. Both
/// produce the same ~200 year time range, ensuring backward pass messages
/// overlap across any reasonable phylogenetic tree time span.
///
/// v0 uses MAX_BRANCH_LENGTH = 4.0 subs/site with non-uniform grids
/// (branch_len_interpolator.py:59). v1 uses uniform grids, so a fixed
/// subs/site limit either wastes resolution for slow rates or is too narrow
/// for fast rates. The time-based limit avoids this.
///
/// At the grid boundary, true Poisson log-likelihood differences are large
/// (~-90 for 23 years at flu rate) but exp(-90) ≈ 3.6e-39 is representable
/// in f64. Values underflow only beyond ~250 years at flu rate, well past
/// this limit.
const MAX_BRANCH_TIME: f64 = 200.0;

fn create_simple_grid(center: f64, one_mutation: f64, n_points: usize, clock_rate: f64) -> Array1<f64> {
  let min_bl = one_mutation * 0.1;
  let peak_max_bl = f64::max(center * 3.0, one_mutation * 10.0);
  let rate_max_bl = MAX_BRANCH_TIME * clock_rate;
  let max_bl = peak_max_bl.max(rate_max_bl);
  Array1::linspace(min_bl, max_bl, n_points)
}
