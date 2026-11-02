use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::timetree::inference::runner::EPS;
use eyre::{Report, WrapErr};
use ndarray::{Array1, ArrayView1, Axis, concatenate, s};
use ndarray_stats::QuantileExt;
use std::iter::successors;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_distribution::DistributionFunction;
use treetime_distribution::NegLog;
use treetime_distribution::rewindow_to_mass;
use treetime_grid::GridFn;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::array::ndarray::{first, last};
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

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

  let PilotGrid {
    grid,
    log_lh,
    right_edge,
  } = pilot_grid(
    contributions,
    indels,
    current_branch_length,
    one_mutation,
    n_grid_points,
    min_bl,
  )?;

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

  let right_boundary = match right_edge {
    RightEdge::Hard => BoundaryBehavior::Hard,
    RightEdge::Soft => {
      let law = SoftTailLaw::fit(&f, Side::Right, DEFAULT_TAIL_FIT_POINTS).wrap_err_with(|| {
        format!("When fitting the branch-length likelihood soft tail near t_max over [{time_min}, {time_max}]")
      })?;
      BoundaryBehavior::Linear(law)
    },
  };

  let distribution_fn = DistributionFunction::from_grid_fn(f)
    .with_left_extrap(left_boundary)?
    .with_right_extrap(right_boundary)?;

  let distribution = Distribution::Function(distribution_fn);
  let distribution = rewindow_to_mass(&distribution, EPS, n_grid_points, max_points)?;
  Ok(Arc::new(distribution))
}

const MIN_BRANCH_LENGTH_MUTATION_FRACTION: f64 = 0.01;
const MAX_BRANCH_LENGTH: f64 = 5.0;
const DECAY_BISECTION_STEPS: usize = 4;
const MAX_PILOT_GROWTH: usize = 10;

fn pilot_grid(
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  center: f64,
  one_mutation: f64,
  n_points: usize,
  min_bl: f64,
) -> Result<PilotGrid, Report> {
  let max_bl = divergence_upper_end(center, one_mutation);
  let pilot = PilotGrid::evaluate(
    Array1::linspace(min_bl, max_bl, n_points),
    RightEdge::Soft,
    contributions,
    indels,
  )?;
  let max_points = n_points * MAX_PILOT_GROWTH;
  match pilot_end(&pilot, max_bl, contributions, indels)? {
    PilotEnd::Keep(right_edge) => Ok(PilotGrid { right_edge, ..pilot }),
    PilotEnd::ExtendTo(max_bl) => pilot.extend_to(max_bl, max_points, contributions, indels),
    PilotEnd::ExtendToCap => {
      let n_points = pilot.points_at_spacing(MAX_BRANCH_LENGTH).min(max_points);
      PilotGrid::evaluate(
        Array1::linspace(min_bl, MAX_BRANCH_LENGTH, n_points),
        RightEdge::Hard,
        contributions,
        indels,
      )
    },
  }
}

fn divergence_upper_end(center: f64, one_mutation: f64) -> f64 {
  const PEAK_BRANCH_LENGTH_MULTIPLE: f64 = 5.0;
  const MIN_GRID_SPAN_MUTATIONS: f64 = 10.0;

  f64::max(
    center * PEAK_BRANCH_LENGTH_MULTIPLE,
    one_mutation * MIN_GRID_SPAN_MUTATIONS,
  )
  .min(MAX_BRANCH_LENGTH)
}

fn pilot_end(
  pilot: &PilotGrid,
  max_bl: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
) -> Result<PilotEnd, Report> {
  let mut probe = DecayProbe {
    contributions,
    indels,
    peak: *pilot.log_lh.max()?,
    floor: *pilot.log_lh.min()?,
    drop: -EPS.ln(),
  };
  if probe.is_decayed(last(&pilot.log_lh)) {
    return Ok(PilotEnd::Keep(RightEdge::Soft));
  }

  let mut below = max_bl;
  for branch_length in outward_probes(max_bl) {
    if probe.decays_at(branch_length)? {
      return Ok(PilotEnd::ExtendTo(probe.bisect(below, branch_length)?));
    }
    below = branch_length;
  }

  Ok(if probe.is_constant() {
    PilotEnd::Keep(RightEdge::Soft)
  } else if max_bl < MAX_BRANCH_LENGTH {
    PilotEnd::ExtendToCap
  } else {
    PilotEnd::Keep(RightEdge::Hard)
  })
}

#[expect(
  variant_size_differences,
  reason = "a pilot end lives only until the grid is chosen, so the size of the larger variant costs nothing"
)]
enum PilotEnd {
  Keep(RightEdge),
  ExtendTo(f64),
  ExtendToCap,
}

fn outward_probes(end: f64) -> impl Iterator<Item = f64> {
  successors(Some(2.0 * end), |x| Some(2.0 * x))
    .take_while(|x| *x < MAX_BRANCH_LENGTH)
    .chain((end < MAX_BRANCH_LENGTH).then_some(MAX_BRANCH_LENGTH))
}

struct PilotGrid {
  grid: Array1<f64>,
  log_lh: Array1<f64>,
  right_edge: RightEdge,
}

impl PilotGrid {
  fn evaluate(
    grid: Array1<f64>,
    right_edge: RightEdge,
    contributions: &[OptimizationContribution],
    indels: IndelTerm,
  ) -> Result<Self, Report> {
    let log_lh = evaluate_log_lh_on_grid(grid.view(), contributions, indels)?;
    Ok(Self {
      grid,
      log_lh,
      right_edge,
    })
  }

  fn spacing(&self) -> f64 {
    (last(&self.grid) - first(&self.grid)) / (self.grid.len() - 1).to_f64_exact()
  }

  fn points_at_spacing(&self, max_bl: f64) -> usize {
    f64_to_usize_saturating(((max_bl - first(&self.grid)) / self.spacing()).ceil()) + 1
  }

  fn extend_to(
    self,
    max_bl: f64,
    max_points: usize,
    contributions: &[OptimizationContribution],
    indels: IndelTerm,
  ) -> Result<Self, Report> {
    let n_points = self.grid.len();
    let min_bl = first(&self.grid);
    let n_extended = self.points_at_spacing(max_bl);
    if n_extended > max_points {
      return Self::evaluate(
        Array1::linspace(min_bl, max_bl, max_points),
        RightEdge::Soft,
        contributions,
        indels,
      );
    }
    let max_bl = min_bl + self.spacing() * (n_extended - 1).to_f64_exact();
    let grid = Array1::linspace(min_bl, max_bl, n_extended);
    let added = evaluate_log_lh_on_grid(grid.slice(s![n_points..]), contributions, indels)?;
    let log_lh = concatenate(Axis(0), &[self.log_lh.view(), added.view()])?;
    Ok(Self {
      grid,
      log_lh,
      right_edge: RightEdge::Soft,
    })
  }
}

#[derive(Clone, Copy, Debug)]
enum RightEdge {
  Soft,
  Hard,
}

struct DecayProbe<'a> {
  contributions: &'a [OptimizationContribution],
  indels: IndelTerm,
  peak: f64,
  floor: f64,
  drop: f64,
}

impl DecayProbe<'_> {
  fn decays_at(&mut self, branch_length: f64) -> Result<bool, Report> {
    let log_lh = evaluate_with_indels_log_lh_only(self.contributions, self.indels, branch_length)?.value();
    self.peak = self.peak.max(log_lh);
    self.floor = self.floor.min(log_lh);
    Ok(self.is_decayed(log_lh))
  }

  fn bisect(&mut self, below: f64, above: f64) -> Result<f64, Report> {
    let mut below = below;
    let mut above = above;
    for _ in 0..DECAY_BISECTION_STEPS {
      let middle = f64::midpoint(below, above);
      if self.decays_at(middle)? {
        above = middle;
      } else {
        below = middle;
      }
    }
    Ok(above)
  }

  fn is_decayed(&self, log_lh: f64) -> bool {
    self.peak - log_lh >= self.drop
  }

  #[expect(
    clippy::float_cmp,
    reason = "a constant likelihood has the same value at every evaluated branch length"
  )]
  fn is_constant(&self) -> bool {
    self.peak == self.floor
  }
}

#[expect(
  clippy::needless_pass_by_value,
  reason = "array views are Copy; clippy misses the Copy impl of ndarray 0.17 views"
)]
fn evaluate_log_lh_on_grid(
  grid: ArrayView1<f64>,
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
