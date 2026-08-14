use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimization_contribution::OptimizationContribution;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;
use treetime_distribution::NegLog;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side};
use treetime_utils::make_report;

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
) -> Result<Arc<Distribution<NegLog>>, Report> {
  debug_assert!(clock_rate > 0.0, "clock_rate must be positive, got {clock_rate:.6e}");
  debug_assert!(gamma > 0.0);

  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points);

  // `create_simple_grid` always returns strictly positive branch lengths
  // (`min_bl = one_mutation * 0.01 > 0`), satisfying the `t > 0` precondition
  // of `poisson_indel_log_lh` when `indel_count > 0`.
  let log_lh: Array1<f64> = grid
    .iter()
    .copied()
    .map(|branch_len| {
      evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, branch_len).map(|log_lh| log_lh.value())
    })
    .collect::<Result<_, _>>()?;
  let max_log_lh = *log_lh.max()?;

  // Store negative-log ordinates directly: `-ln(p / p_peak) = ln(p_peak) - ln(p) = max_log_lh -
  // log_lh`. The peak ordinate is `0`, matching the `NegLog` peak-normalization convention, and no
  // `exp` round-trip is taken.
  let neg_log = log_lh.mapv(|value| max_log_lh - value);

  // Convert branch length grid to time grid: time = branch_length / (clock_rate * gamma)
  // gamma > 1 means faster evolution, so same substitutions correspond to shorter time
  let effective_clock_rate = clock_rate * gamma;
  let time_min = grid[0] / effective_clock_rate;
  let time_max = grid[grid.len() - 1] / effective_clock_rate;

  let distribution_fn = DistributionFunction::from_range_values((time_min, time_max), neg_log)?;

  // Fit the hard-boundary approach law on the left side (near t=0). The grid stores neg-log
  // ordinates and starts at min_bl / clock_rate > 0. The fit detects the regime from the data:
  // for n >= 1 mutations the neg-log density diverges logarithmically (power-law approach),
  // for n = 0 mutations it is exactly linear (the mode sits on the boundary).
  //
  // A grid too degenerate to fit the law (fewer than two finite innermost points) is an error, not
  // a silent flat boundary: the gap `[0, t_first)` carries real density (maximal at t=0 for a
  // zero-mutation branch), so a missing law would zero it and lose that mass -- the class of defect
  // that motivated the hard-approach law.
  let approach =
    HardApproachLaw::fit(distribution_fn.grid_fn(), 0.0, Side::Left, DEFAULT_TAIL_FIT_POINTS).ok_or_else(|| {
      make_report!(
        "Branch-length likelihood grid over [{time_min}, {time_max}] is too degenerate to fit a \
         hard-boundary approach law near t=0 (fewer than two finite innermost points)"
      )
    })?;
  let distribution_fn = distribution_fn.with_left_extrap(BoundaryBehavior::HardApproach(approach))?;

  Ok(Arc::new(Distribution::Function(distribution_fn)))
}

fn create_simple_grid(center: f64, one_mutation: f64, n_points: usize) -> Array1<f64> {
  // Grid floor as a fraction of one mutation's branch length. Keeps `t > 0` for
  // the Poisson indel term while sitting well below any resolvable branch length.
  const MIN_BRANCH_LENGTH_MUTATION_FRACTION: f64 = 0.01;
  // Grid upper bound as a multiple of the current ML branch length, so grid
  // resolution concentrates around the likelihood peak.
  const PEAK_BRANCH_LENGTH_MULTIPLE: f64 = 5.0;
  // Minimum grid span as a multiple of one mutation, for near-zero branch lengths.
  const MIN_GRID_SPAN_MUTATIONS: f64 = 10.0;
  // Absolute cap on the grid upper bound (substitutions/site).
  const MAX_BRANCH_LENGTH: f64 = 5.0;

  let min_bl = one_mutation * MIN_BRANCH_LENGTH_MUTATION_FRACTION;
  let peak_max_bl = f64::max(
    center * PEAK_BRANCH_LENGTH_MULTIPLE,
    one_mutation * MIN_GRID_SPAN_MUTATIONS,
  );
  let max_bl = peak_max_bl.min(MAX_BRANCH_LENGTH);
  Array1::linspace(min_bl, max_bl, n_points)
}
