use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimization_contribution::OptimizationContribution;
use crate::timetree::inference::runner::EPS;
use eyre::{Report, WrapErr};
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;
use treetime_distribution::NegLog;
use treetime_distribution::rewindow_to_mass;
use treetime_grid::GridFn;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side, SoftTailLaw};
use treetime_utils::array::ndarray::{first, last};

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

  let log_lh = evaluate_log_lh_on_grid(&grid, contributions, indel_count, indel_rate)?;

  let (log_lh, max_log_lh) = peak_normalize_neg_log(&log_lh)?;

  let TimeRange { time_min, time_max } = branch_length_grid_to_time_range(&grid, clock_rate, gamma);

  let f = GridFn::from_range_values((time_min, time_max), log_lh)?;

  let y_hard = branch_length_boundary_ordinate(contributions, indel_count, indel_rate, max_log_lh)?;
  // A finite boundary ordinate selects the straight-line fit; its absence (a forced substitution or
  // an indel, whose neg-log diverges at t=0) selects the power-law fit.
  let left_boundary = match y_hard {
    Some(y_hard) => HardApproachLaw::fit_linear(&f, 0.0, Side::Left, y_hard),
    None => HardApproachLaw::fit_log_power_law(&f, 0.0, Side::Left, DEFAULT_TAIL_FIT_POINTS),
  }
  .wrap_err_with(|| {
    format!("When building the branch-length likelihood hard boundary near t=0 over [{time_min}, {time_max}]")
  })?;

  let right_boundary = SoftTailLaw::fit(&f, Side::Right, DEFAULT_TAIL_FIT_POINTS).wrap_err_with(|| {
    format!("When fitting the branch-length likelihood soft tail near t_max over [{time_min}, {time_max}]")
  })?;

  let distribution_fn = DistributionFunction::from_grid_fn(f)
    .with_left_extrap(BoundaryBehavior::HardApproach(left_boundary))?
    .with_right_extrap(BoundaryBehavior::Linear(right_boundary))?;

  // The heuristic pilot grid above is tail-complete (left HardApproach, right Linear). Re-window it by
  // probability mass so the stored grid holds >= 1 - 2*EPS of the evaluated likelihood, sized by mass
  // rather than by heuristic multiples of the branch length (design D1). A zero-mutation branch keeps
  // its lower edge on the finite hard bound at t = 0.
  let distribution = Distribution::Function(distribution_fn);
  let distribution = rewindow_to_mass(&distribution, EPS, n_grid_points)?;
  Ok(Arc::new(distribution))
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

/// Evaluate the combined substitution + Poisson-indel log-likelihood at each grid branch length.
///
/// `create_simple_grid` always returns strictly positive branch lengths
/// (`min_bl = one_mutation * 0.01 > 0`), satisfying the `t > 0` precondition
/// of `poisson_indel_log_lh` when `indel_count > 0`.
fn evaluate_log_lh_on_grid(
  grid: &Array1<f64>,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
) -> Result<Array1<f64>, Report> {
  grid
    .iter()
    .copied()
    .map(|branch_len| {
      evaluate_with_indels_log_lh_only(contributions, indel_count, indel_rate, branch_len).map(|lh| lh.value())
    })
    .collect()
}

/// Peak-normalize log-likelihood ordinates into the `NegLog` convention.
///
/// Returns the neg-log ordinates `max_log_lh - log_lh` (peak ordinate `0`, all others `>= 0`)
/// alongside the peak log-likelihood `max_log_lh`, so callers can normalize additional off-grid
/// ordinates (such as the boundary at branch length 0) against the same peak. The subtraction
/// stores `-ln(p / p_peak) = ln(p_peak) - ln(p)` directly, taking no `exp` round-trip.
fn peak_normalize_neg_log(log_lh: &Array1<f64>) -> Result<(Array1<f64>, f64), Report> {
  let max_log_lh = *log_lh.max()?;
  let neg_log = log_lh.mapv(|value| max_log_lh - value);
  Ok((neg_log, max_log_lh))
}

/// Convert a branch-length grid (substitutions/site) to the corresponding time range.
///
/// `time = branch_length / (clock_rate * gamma)`. A `gamma > 1` (faster evolution) maps the same
/// substitution count to a shorter time, compressing the range.
fn branch_length_grid_to_time_range(grid: &Array1<f64>, clock_rate: f64, gamma: f64) -> TimeRange {
  let effective_clock_rate = clock_rate * gamma;
  TimeRange {
    time_min: first(grid) / effective_clock_rate,
    time_max: last(grid) / effective_clock_rate,
  }
}

struct TimeRange {
  time_min: f64,
  time_max: f64,
}

/// Boundary neg-log ordinate at `t = 0`, peak-normalized against `max_log_lh`.
///
/// Classifies the hard boundary of the branch-length density: `Some(max_log_lh - log_lh(0))` when
/// the density is finite there (a zero-mutation branch, whose neg-log approaches the boundary as a
/// straight line), and `None` when it diverges (a forced substitution or an indel, whose neg-log
/// approaches as a power law). The `Some` case supplies the boundary ordinate to
/// `HardApproachLaw::fit_linear`; the `None` case selects `HardApproachLaw::fit_log_power_law`.
///
/// An indel branch cannot be evaluated at `t = 0` (`k*ln(0) = -inf`, and `poisson_indel_log_lh`
/// rejects `t = 0` for `k > 0`), so `indel_count > 0` diverges without evaluation.
fn branch_length_boundary_ordinate(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  max_log_lh: f64,
) -> Result<Option<f64>, Report> {
  if indel_count > 0 {
    return Ok(None);
  }
  let boundary_log_lh = evaluate_with_indels_log_lh_only(contributions, 0, indel_rate, 0.0)?.value();
  Ok(boundary_log_lh.is_finite().then_some(max_log_lh - boundary_log_lh))
}
