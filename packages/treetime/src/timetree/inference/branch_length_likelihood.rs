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
  indels: IndelTerm,
  current_branch_length: f64,
  one_mutation: f64,
  n_grid_points: usize,
  branch_rate: f64,
  max_points: MaxGridPoints,
) -> Result<Arc<Distribution<NegLog>>, Report> {
  debug_assert!(branch_rate > 0.0, "branch_rate must be positive, got {branch_rate:.6e}");

  let finite_boundary = branch_length_has_finite_boundary(contributions, indels)?;
  let min_bl = if finite_boundary {
    0.0
  } else {
    one_mutation * MIN_BRANCH_LENGTH_MUTATION_FRACTION
  };

  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points, min_bl);

  let log_lh = evaluate_log_lh_on_grid(&grid, contributions, indels)?;

  let log_lh = peak_normalize_neg_log(&log_lh)?;

  let TimeRange { time_min, time_max } = branch_length_grid_to_time_range(&grid, branch_rate);

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
  let distribution = rewindow_to_mass(&distribution, EPS, n_grid_points, max_points)?;
  Ok(Arc::new(distribution))
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

fn branch_length_grid_to_time_range(grid: &Array1<f64>, branch_rate: f64) -> TimeRange {
  TimeRange {
    time_min: first(grid) / branch_rate,
    time_max: last(grid) / branch_rate,
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
