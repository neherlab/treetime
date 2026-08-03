use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimization_contribution::OptimizationContribution;
use crate::timetree::inference::runner::BRANCH_GRID_SIZE;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;
use treetime_distribution::TailSide;

/// Probability floor, relative to the peak, below which the extrapolated long-branch tail is
/// numerically negligible and extension stops.
const TAIL_REL_FLOOR: f64 = 1e-10;
/// Maximum multiple of the base grid size to which the long-branch tail extension may grow the
/// branch-length grid. The relative floor terminates earlier for informative branches, so this
/// only bounds convolution cost and near-flat (weakly informative) tails. The 10x budget matches
/// the grid extent whose corrected node dates were validated in the sensitivity study behind
/// kb/decisions/distribution-tails-and-arithmetic.md.
const TAIL_MAX_GRID_GROWTH: usize = 10;
/// Default grid upper bound as a multiple of the current ML branch length, so grid resolution
/// concentrates around the likelihood peak.
const PEAK_BRANCH_LENGTH_MULTIPLE: f64 = 5.0;
/// Default absolute cap on the grid upper bound (substitutions/site).
const MAX_BRANCH_LENGTH: f64 = 5.0;

/// Runtime configuration for the branch-length likelihood grid and its long-branch tail extension.
///
/// The knobs govern the accuracy/cost tradeoff of the branch-length distribution used in timetree
/// inference: `grid_size` is the peak resolution; `peak_extent_multiple` and `max_branch_length`
/// set how far the computed (real-likelihood) grid reaches before truncation; and
/// `tail_max_grid_growth` with `tail_rel_floor` govern how far the extrapolated log-linear far-tail
/// extension reaches beyond that. `Default` reproduces the values validated in
/// `kb/decisions/distribution-tails-and-arithmetic.md`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct BranchGridConfig {
  /// Number of uniform grid points across the peak-scaled branch-length interval.
  pub grid_size: usize,
  /// Grid upper bound as a multiple of the ML branch length (peak-scaled real-likelihood extent).
  pub peak_extent_multiple: f64,
  /// Absolute cap on the grid upper bound (substitutions/site).
  pub max_branch_length: f64,
  /// Maximum multiple of `grid_size` to which the extrapolated long-branch tail may grow the grid.
  pub tail_max_grid_growth: usize,
  /// Probability floor, relative to the peak, below which tail extension stops.
  pub tail_rel_floor: f64,
}

impl Default for BranchGridConfig {
  fn default() -> Self {
    Self {
      grid_size: BRANCH_GRID_SIZE,
      peak_extent_multiple: PEAK_BRANCH_LENGTH_MULTIPLE,
      max_branch_length: MAX_BRANCH_LENGTH,
      tail_max_grid_growth: TAIL_MAX_GRID_GROWTH,
      tail_rel_floor: TAIL_REL_FLOOR,
    }
  }
}

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
  grid_config: &BranchGridConfig,
  clock_rate: f64,
  gamma: f64,
) -> Result<Arc<Distribution>, Report> {
  debug_assert!(clock_rate > 0.0, "clock_rate must be positive, got {clock_rate:.6e}");
  debug_assert!(gamma > 0.0);

  let grid = create_simple_grid(
    current_branch_length,
    one_mutation,
    grid_config.grid_size,
    grid_config.peak_extent_multiple,
    grid_config.max_branch_length,
  );

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
  let max_log_lh = log_lh.max()?;

  let normalized_prob = (&log_lh - *max_log_lh).exp();

  // `create_simple_grid` hard-truncates the long-branch side, but the branch-length likelihood
  // decays exponentially for long branches. Continue that decay as a log-linear tail so the
  // far-past reach of the backward message reflects the physical decay instead of a flat Constant
  // floor, which otherwise biases deep internal-node dates older
  // (kb/decisions/distribution-tails-and-arithmetic.md).
  let branch_dist = DistributionFunction::from_arrays(&grid, normalized_prob)?.extend_log_linear_tail(
    TailSide::Right,
    grid_config.tail_rel_floor,
    (grid_config.tail_max_grid_growth - 1) * grid_config.grid_size,
  )?;

  // Convert branch length grid to time grid: time = branch_length / (clock_rate * gamma)
  // gamma > 1 means faster evolution, so same substitutions correspond to shorter time. The
  // conversion is a uniform axis rescale, so it preserves the extended grid's uniform spacing.
  let effective_clock_rate = clock_rate * gamma;
  let time_min = branch_dist.x_min() / effective_clock_rate;
  let time_dx = branch_dist.dx() / effective_clock_rate;

  let distribution_fn = DistributionFunction::from_start_dx_values(time_min, time_dx, branch_dist.y().clone())?;
  Ok(Arc::new(Distribution::Function(distribution_fn)))
}

fn create_simple_grid(
  center: f64,
  one_mutation: f64,
  n_points: usize,
  peak_extent_multiple: f64,
  max_branch_length: f64,
) -> Array1<f64> {
  // Grid floor as a fraction of one mutation's branch length. Keeps `t > 0` for
  // the Poisson indel term while sitting well below any resolvable branch length.
  const MIN_BRANCH_LENGTH_MUTATION_FRACTION: f64 = 0.01;
  // Minimum grid span as a multiple of one mutation, for near-zero branch lengths.
  const MIN_GRID_SPAN_MUTATIONS: f64 = 10.0;

  let min_bl = one_mutation * MIN_BRANCH_LENGTH_MUTATION_FRACTION;
  let peak_max_bl = f64::max(center * peak_extent_multiple, one_mutation * MIN_GRID_SPAN_MUTATIONS);
  let max_bl = peak_max_bl.min(max_branch_length);
  Array1::linspace(min_bl, max_bl, n_points)
}
