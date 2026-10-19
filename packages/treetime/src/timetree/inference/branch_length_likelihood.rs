use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimize::contribution::OptimizationContribution;
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
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::array::ndarray::{first, last};

pub(crate) fn compute_branch_length_distribution(
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  current_branch_length: f64,
  one_mutation: f64,
  n_grid_points: usize,
  clock_rate: f64,
  gamma: f64,
  max_points: MaxGridPoints,
) -> Result<Arc<Distribution<NegLog>>, Report> {
  debug_assert!(clock_rate > 0.0, "clock_rate must be positive, got {clock_rate:.6e}");
  debug_assert!(gamma > 0.0);

  let indels = IndelTerm::new(indel_count, indel_rate);
  let finite_boundary = branch_length_has_finite_boundary(contributions, indels)?;
  let min_bl = if finite_boundary {
    0.0
  } else {
    one_mutation * MIN_BRANCH_LENGTH_MUTATION_FRACTION
  };

  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points, min_bl);

  let started = std::time::Instant::now();
  let log_lh = evaluate_log_lh_on_grid(&grid, contributions, indels)?;
  experiment_timers::add(1, started.elapsed());
  let started = std::time::Instant::now();

  let log_lh = peak_normalize_neg_log(&log_lh)?;

  let TimeRange { time_min, time_max } = branch_length_grid_to_time_range(&grid, clock_rate, gamma);

  let f = GridFn::from_range_values((time_min, time_max), log_lh)?;

  let left_boundary = if finite_boundary {
    BoundaryBehavior::Hard
  } else {
    let law = HardApproachLaw::fit(&f, 0.0, Side::Left, DEFAULT_TAIL_FIT_POINTS).wrap_err_with(|| {
      format!("When building the branch-length likelihood hard boundary near t=0 over [{time_min}, {time_max}]")
    })?;
    BoundaryBehavior::HardApproach(law)
  };

  let right_boundary = SoftTailLaw::fit(&f, Side::Right, DEFAULT_TAIL_FIT_POINTS).wrap_err_with(|| {
    format!("When fitting the branch-length likelihood soft tail near t_max over [{time_min}, {time_max}]")
  })?;

  let distribution_fn = DistributionFunction::from_grid_fn(f)
    .with_left_extrap(left_boundary)?
    .with_right_extrap(BoundaryBehavior::Linear(right_boundary))?;

  let distribution = Distribution::Function(distribution_fn);
  experiment_timers::add(2, started.elapsed());
  let started = std::time::Instant::now();
  let distribution = rewindow_to_mass(&distribution, EPS, n_grid_points, max_points)?;
  experiment_timers::add(3, started.elapsed());
  Ok(Arc::new(distribution))
}

/// EXPERIMENT: time split of the branch distributions, printed per round with `TT_GRID_STATS`
pub(crate) mod experiment_timers {
  use std::sync::atomic::{AtomicU64, Ordering};
  use std::time::Duration;

  static NANOS: [AtomicU64; 4] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];

  pub(crate) fn add(slot: usize, elapsed: Duration) {
    NANOS[slot].fetch_add(u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX), Ordering::Relaxed);
  }

  pub(crate) fn report() {
    if std::env::var_os("TT_GRID_STATS").is_none() {
      return;
    }
    let secs: Vec<f64> = NANOS.iter().map(|n| n.swap(0, Ordering::Relaxed) as f64 * 1e-9).collect();
    eprintln!(
      "BRANCH_STATS contributions={:.3} grid_eval={:.3} tails={:.3} rewindow={:.3}",
      secs[0], secs[1], secs[2], secs[3]
    );
  }
}

const MIN_BRANCH_LENGTH_MUTATION_FRACTION: f64 = 0.01;

fn create_simple_grid(center: f64, one_mutation: f64, n_points: usize, min_bl: f64) -> Array1<f64> {
  const PEAK_BRANCH_LENGTH_MULTIPLE: f64 = 5.0;
  const MIN_GRID_SPAN_MUTATIONS: f64 = 10.0;
  const MAX_BRANCH_LENGTH: f64 = 5.0;

  let peak_max_bl = f64::max(
    center * PEAK_BRANCH_LENGTH_MULTIPLE,
    one_mutation * MIN_GRID_SPAN_MUTATIONS,
  );
  let max_bl = peak_max_bl.min(MAX_BRANCH_LENGTH);
  Array1::linspace(min_bl, max_bl, n_points)
}

fn evaluate_log_lh_on_grid(
  grid: &Array1<f64>,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
) -> Result<Array1<f64>, Report> {
  grid
    .iter()
    .copied()
    .map(|branch_len| evaluate_with_indels_log_lh_only(contributions, indels, branch_len).map(|lh| lh.value()))
    .collect()
}

fn peak_normalize_neg_log(log_lh: &Array1<f64>) -> Result<Array1<f64>, Report> {
  let max_log_lh = *log_lh.max()?;
  Ok(log_lh.mapv(|value| max_log_lh - value))
}

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

fn branch_length_has_finite_boundary(
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
) -> Result<bool, Report> {
  if indels.requires_positive_length() {
    return Ok(false);
  }
  let boundary_log_lh = evaluate_with_indels_log_lh_only(contributions, indels, 0.0)?.value();
  Ok(boundary_log_lh.is_finite())
}
