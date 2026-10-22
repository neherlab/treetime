use argmin::core::{CostFunction, Error, Executor, TerminationReason, TerminationStatus};
use argmin::solver::brent::BrentRoot;
use eyre::Report;
use itertools::Itertools;
use ndarray::{Array1, ArrayView1, s};
use ordered_float::OrderedFloat;
use std::iter::once;
use treetime_distribution::{
  BoundaryBehavior, Distribution, DistributionFunction, MassFractions, MassProfile, NegLog, TailMass,
};
use treetime_grid::{EdgeCell, GridEdge, HardApproachLaw, HardCell, Side};
use treetime_utils::array::ndarray::partition_point;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error, make_report};

const LEVEL_BRACKET_DOUBLINGS: usize = 64;
const LEVEL_SOLVE_MAX_ITERS: u64 = 200;
const REGION_MASS_ROUNDING_ULPS: f64 = 64.0;

pub fn max_posterior_region(
  distribution: &Distribution<NegLog>,
  fraction: f64,
) -> Result<Option<MaxPosteriorRegion>, Report> {
  if !(fraction > 0.0 && fraction < 1.0) {
    return make_error!("Posterior mass fraction must be in (0, 1), got {fraction}");
  }
  let region = match Posterior::of(distribution)? {
    Posterior::Empty => return Ok(None),
    Posterior::Point(_) => return Ok(Some(MaxPosteriorRegion::Point)),
    Posterior::Range(start, end) => PosteriorRegion {
      lower: start + (1.0 - fraction) * (end - start),
      upper: end,
      peak: PeakLocation::Latest,
    },
    Posterior::Function(function) => function_max_posterior_region(&function, fraction)?,
  };
  Ok(Some(MaxPosteriorRegion::Grid(region)))
}

pub(crate) fn posterior_median(distribution: &Distribution<NegLog>) -> Result<Option<f64>, Report> {
  match Posterior::of(distribution)? {
    Posterior::Empty => Ok(None),
    Posterior::Point(t) => Ok(Some(t)),
    Posterior::Range(start, end) => Ok(Some(f64::midpoint(start, end))),
    Posterior::Function(function) => {
      let profile = MassProfile::new(&function)?;
      Ok(Some(posterior_fractions(&profile)?.quantile(0.5)))
    },
  }
}

pub(crate) fn posterior_survival(
  distribution: &Distribution<NegLog>,
  levels: &[SurvivalLevel],
) -> Result<Option<PosteriorSurvival>, Report> {
  match Posterior::of(distribution)? {
    Posterior::Empty => Ok(None),
    Posterior::Point(t) => Ok(Some(PosteriorSurvival::point(t))),
    Posterior::Range(start, end) => Ok(Some(PosteriorSurvival {
      t: vec![start, end],
      survival: vec![1.0, 0.0],
    })),
    Posterior::Function(function) => function_survival(&function, levels).map(Some),
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum SurvivalLevel {
  Below(f64),
  Above(f64),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaxPosteriorRegion {
  Point,
  Grid(PosteriorRegion),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PosteriorRegion {
  pub lower: f64,
  pub upper: f64,
  pub peak: PeakLocation,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeakLocation {
  Earliest,
  Interior,
  Latest,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PosteriorSurvival {
  pub t: Vec<f64>,
  pub survival: Vec<f64>,
}

impl PosteriorSurvival {
  fn point(t: f64) -> Self {
    Self {
      t: vec![t],
      survival: vec![0.0],
    }
  }

  pub(crate) fn push(&mut self, t: f64, survival: f64) {
    let later = self.t.last().is_none_or(|&last| t > last);
    let lower = self.survival.last().is_none_or(|&last| survival < last);
    if later && lower {
      self.t.push(t);
      self.survival.push(survival);
    }
  }

  pub(crate) fn with_unit_ends(mut self) -> Self {
    if let [first, .., last] = self.survival.as_mut_slice() {
      *first = 1.0;
      *last = 0.0;
    } else {
      self.survival = vec![0.0];
    }
    self
  }
}

enum Posterior {
  Empty,
  Point(f64),
  Range(f64, f64),
  Function(DistributionFunction<f64, NegLog>),
}

impl Posterior {
  fn of(distribution: &Distribution<NegLog>) -> Result<Self, Report> {
    let normalized = distribution.normalize()?;
    let bounds = normalized.time_bounds();
    Ok(match (normalized, bounds) {
      (Distribution::Function(function), _) => Self::Function(function),
      (Distribution::Point(_), Some((t, _))) => Self::Point(t),
      (Distribution::Range(_), Some((start, end))) => Self::Range(start, end),
      (Distribution::Empty, _) => Self::Empty,
      (Distribution::Formula(_), _) => {
        return make_internal_error!("Normalization left a formula distribution undiscretized");
      },
      (Distribution::Point(_) | Distribution::Range(_), None) => {
        return make_internal_error!("A point or range distribution has no time bounds");
      },
    })
  }
}

fn posterior_fractions(profile: &MassProfile) -> Result<MassFractions<'_>, Report> {
  match profile.fractions(TailMass::ZeroWhereUndefined) {
    Some(fractions) => Ok(fractions),
    None => make_error!(
      "Posterior on [{}, {}] has no probability mass",
      profile.x_min(),
      profile.x_max()
    ),
  }
}

fn function_survival(
  function: &DistributionFunction<f64, NegLog>,
  levels: &[SurvivalLevel],
) -> Result<PosteriorSurvival, Report> {
  let profile = MassProfile::new(function)?;
  let fractions = posterior_fractions(&profile)?;
  let time_at = |level: SurvivalLevel| match level {
    SurvivalLevel::Below(mass) => fractions.quantile(mass),
    SurvivalLevel::Above(mass) => fractions.survival_inverse(mass),
  };
  for &level in levels.first().into_iter().chain(levels.last()) {
    let t = time_at(level);
    if t.is_finite() {
      continue;
    }
    let (SurvivalLevel::Below(mass) | SurvivalLevel::Above(mass)) = level;
    let (side, tail) = if t > 0.0 {
      (Side::Right, function.right_extrap())
    } else {
      (Side::Left, function.left_extrap())
    };
    return match tail {
      BoundaryBehavior::Linear(law) => make_error!(
        "Posterior on [{}, {}] has a tail of slope {:e} on the {side:?} side, which holds mass {mass:e} beyond the largest finite time",
        function.x_min(),
        function.x_max(),
        law.slope
      ),
      BoundaryBehavior::Hard
      | BoundaryBehavior::HardApproach(_)
      | BoundaryBehavior::HardEdgeCell(_)
      | BoundaryBehavior::Error => {
        make_internal_error!("Posterior quantile {mass:e} on the {side:?} side without a soft tail is not finite")
      },
    };
  }
  let mut survival = PosteriorSurvival {
    t: Vec::with_capacity(levels.len()),
    survival: Vec::with_capacity(levels.len()),
  };
  for &level in levels {
    let t = time_at(level);
    survival.push(t, fractions.survival(t));
  }
  Ok(survival.with_unit_ends())
}

fn function_max_posterior_region(
  function: &DistributionFunction<f64, NegLog>,
  fraction: f64,
) -> Result<PosteriorRegion, Report> {
  let peak = latest_peak_index(function.y().view());
  let origin = function.x_min() + function.dx() * peak.to_f64_exact();
  let local = translated(function, -origin)?;
  let profile = MassProfile::new(&local)?;
  let Some(fractions) = profile.fractions(TailMass::ZeroWhereUndefined) else {
    return make_error!(
      "Posterior on [{}, {}] has no probability mass",
      function.x_min(),
      function.x_max()
    );
  };
  let local_region = local_max_posterior_region(&local, &fractions, peak, fraction)?;
  let region = PosteriorRegion {
    lower: origin + local_region.lower,
    upper: origin + local_region.upper,
    peak: local_region.peak,
  };
  check_region_mass(&fractions, origin, region, fraction)?;
  Ok(region)
}

fn translated(
  function: &DistributionFunction<f64, NegLog>,
  shift: f64,
) -> Result<DistributionFunction<f64, NegLog>, Report> {
  let tail = |behavior: BoundaryBehavior| match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(HardApproachLaw {
      t_hard: law.t_hard + shift,
      ..law
    }),
    BoundaryBehavior::HardEdgeCell(law) => BoundaryBehavior::HardEdgeCell(law.shift_arg(shift)),
    other => other,
  };
  DistributionFunction::from_start_dx_values(function.x_min() + shift, function.dx(), function.y().clone())?
    .with_left_extrap(tail(function.left_extrap()))?
    .with_right_extrap(tail(function.right_extrap()))
}

fn local_max_posterior_region(
  function: &DistributionFunction<f64, NegLog>,
  fractions: &MassFractions<'_>,
  peak: usize,
  fraction: f64,
) -> Result<PosteriorRegion, Report> {
  let n_points = function.y().len();
  let (x_min, x_max) = (function.x_min(), function.x_max());
  let carries_mass = BoundaryMass {
    left: fractions.cdf(x_min) > 0.0,
    right: fractions.survival(x_max) > 0.0,
  };
  let cells = EdgeCells::of(function)?;
  let peak_cell = cells.peak_cell(function.y()[peak]);
  if peak_cell.is_none() && peak == n_points - 1 && !carries_mass.right {
    return Ok(PosteriorRegion {
      lower: fractions.survival_inverse(fraction),
      upper: x_max,
      peak: PeakLocation::Latest,
    });
  }
  if peak_cell.is_none() && peak == 0 && !carries_mass.left {
    return Ok(PosteriorRegion {
      lower: x_min,
      upper: fractions.quantile(fraction),
      peak: PeakLocation::Earliest,
    });
  }

  let (earlier, later) = match peak_cell {
    None => grid_peak_sides(function, peak, carries_mass, cells)?,
    Some((side, cell, mode)) => cell_peak_sides(function, side, cell, mode, carries_mass, cells)?,
  };
  let region = LevelSolve {
    fractions,
    earlier,
    later,
    fraction,
  }
  .solve()?;
  Ok(PosteriorRegion {
    lower: region.lower,
    upper: region.upper,
    peak: PeakLocation::Interior,
  })
}

#[derive(Clone, Copy)]
struct BoundaryMass {
  left: bool,
  right: bool,
}

#[derive(Clone, Copy)]
struct EdgeCells {
  left: Option<EdgeCell>,
  right: Option<EdgeCell>,
}

impl EdgeCells {
  fn of(function: &DistributionFunction<f64, NegLog>) -> Result<Self, Report> {
    let last_cell = function.y().len() - 2;
    let left = function.edge_cell(0)?.filter(|cell| cell.edge.t <= cell.anchor.t);
    let right = function
      .edge_cell(last_cell)?
      .filter(|cell| cell.edge.t > cell.anchor.t);
    Ok(Self { left, right })
  }

  fn peak_cell(self, grid_peak: f64) -> Option<(Side, EdgeCell, GridEdge)> {
    [(Side::Left, self.left), (Side::Right, self.right)]
      .into_iter()
      .filter_map(|(side, cell)| cell.and_then(|cell| cell.mode().map(|mode| (side, cell, mode))))
      .filter(|(_, _, mode)| mode.y < grid_peak)
      .min_by_key(|(_, _, mode)| OrderedFloat(mode.y))
  }
}

fn grid_peak_sides(
  function: &DistributionFunction<f64, NegLog>,
  peak: usize,
  carries_mass: BoundaryMass,
  cells: EdgeCells,
) -> Result<(LevelSide, LevelSide), Report> {
  let (t, y) = (grid_t(function), function.y());
  let n_points = y.len();
  let earlier_segments = chain(
    cells.left.filter(|_| peak > 0).map(|cell| (peak, cell)),
    cells.right.filter(|_| peak == n_points - 1).map(|cell| (1, cell)),
  );
  let later_segments = chain(
    cells
      .right
      .filter(|_| peak + 1 < n_points)
      .map(|cell| (n_points - 1 - peak, cell)),
    cells.left.filter(|_| peak == 0).map(|cell| (1, cell)),
  );
  let earlier = LevelSide::new(
    t.slice(s![..=peak;-1]),
    y.slice(s![..=peak;-1]),
    Side::Left,
    carries_mass.left.then_some(function.left_extrap()),
    first_hard_cell(function, (0..peak).rev(), |cell| cell)?,
    earlier_segments.collect(),
  );
  let later = LevelSide::new(
    t.slice(s![peak..]),
    y.slice(s![peak..]),
    Side::Right,
    carries_mass.right.then_some(function.right_extrap()),
    first_hard_cell(function, peak..n_points - 1, |cell| cell + 1)?,
    later_segments.collect(),
  );
  Ok((earlier, later))
}

fn cell_peak_sides(
  function: &DistributionFunction<f64, NegLog>,
  side: Side,
  cell: EdgeCell,
  mode: GridEdge,
  carries_mass: BoundaryMass,
  cells: EdgeCells,
) -> Result<(LevelSide, LevelSide), Report> {
  let (t, y) = (grid_t(function), function.y());
  let n_points = y.len();
  let with_mode = |values: ArrayView1<f64>, at: f64| chain(once(at), values.iter().copied()).collect::<Array1<f64>>();
  let (edge_index, inward, inward_hard) = match side {
    Side::Left => (0, s![1..], first_hard_cell(function, 1..n_points - 1, |cell| cell + 1)?),
    Side::Right => (
      n_points - 1,
      s![..n_points - 1;-1],
      first_hard_cell(function, (0..n_points - 2).rev(), |cell| cell)?,
    ),
  };
  let other = match side {
    Side::Left => cells.right,
    Side::Right => cells.left,
  };
  let outward = LevelSide::new(
    with_mode(t.slice(s![edge_index..=edge_index]), mode.t).view(),
    with_mode(y.slice(s![edge_index..=edge_index]), mode.y).view(),
    side,
    None,
    None,
    vec![(1, cell)],
  );
  let inward_side = match side {
    Side::Left => Side::Right,
    Side::Right => Side::Left,
  };
  let inward_tail = match side {
    Side::Left => carries_mass.right.then_some(function.right_extrap()),
    Side::Right => carries_mass.left.then_some(function.left_extrap()),
  };
  let inward = LevelSide::new(
    with_mode(t.slice(inward), mode.t).view(),
    with_mode(y.slice(inward), mode.y).view(),
    inward_side,
    inward_tail,
    inward_hard,
    chain(
      once((1, cell)),
      other.filter(|_| n_points > 2).map(|other| (n_points - 1, other)),
    )
    .collect(),
  );
  Ok(match side {
    Side::Left => (outward, inward),
    Side::Right => (inward, outward),
  })
}

fn grid_t(function: &DistributionFunction<f64, NegLog>) -> Array1<f64> {
  Array1::from_shape_fn(function.y().len(), |idx| {
    function.x_min() + function.dx() * idx.to_f64_exact()
  })
}

fn first_hard_cell(
  function: &DistributionFunction<f64, NegLog>,
  cells: impl Iterator<Item = usize>,
  outer_point: impl Fn(usize) -> usize,
) -> Result<Option<HardCell>, Report> {
  let y = function.y();
  match cells.into_iter().find(|&cell| y[outer_point(cell)].is_infinite()) {
    Some(cell) => function.hard_cell(cell),
    None => Ok(None),
  }
}

fn check_region_mass(
  fractions: &MassFractions<'_>,
  origin: f64,
  region: PosteriorRegion,
  fraction: f64,
) -> Result<(), Report> {
  let total = fractions.total();
  let (lower, upper) = (region.lower - origin, region.upper - origin);
  let mass = fractions.mass_between(lower, upper)? / total;
  let rounding_mass = |local: f64, absolute: f64| -> Result<f64, Report> {
    let half_ulp = ulp(absolute) / 2.0;
    Ok(fractions.mass_between(local - half_ulp, local + half_ulp)? / total)
  };
  let bound = rounding_mass(lower, region.lower)?
    + rounding_mass(upper, region.upper)?
    + REGION_MASS_ROUNDING_ULPS * f64::EPSILON;
  if (mass - fraction).abs() > bound {
    return make_internal_error!(
      "Maximum posterior region [{}, {}] holds the fraction {mass} instead of {fraction}, beyond the rounding bound {bound}",
      region.lower,
      region.upper
    );
  }
  Ok(())
}

fn ulp(x: f64) -> f64 {
  if x.is_finite() {
    x.abs().next_up() - x.abs()
  } else {
    0.0
  }
}

#[expect(
  clippy::needless_pass_by_value,
  reason = "array views are Copy; clippy misses the Copy impl of ndarray 0.17 views"
)]
fn latest_peak_index(neglog: ArrayView1<f64>) -> usize {
  neglog
    .iter()
    .enumerate()
    .rev()
    .min_by_key(|(_, y)| OrderedFloat(**y))
    .map_or(0, |(i, _)| i)
}

struct LevelSolve<'a> {
  fractions: &'a MassFractions<'a>,
  earlier: LevelSide,
  later: LevelSide,
  fraction: f64,
}

impl LevelSolve<'_> {
  fn solve(&self) -> Result<Interval, Report> {
    let levels = self
      .earlier
      .envelope
      .iter()
      .merge(&self.later.envelope)
      .copied()
      .filter(|level| level.is_finite())
      .dedup()
      .chain(once(f64::INFINITY));

    let mut below = 0.0;
    for level in levels {
      let under = self.region(level, LevelBound::Below);
      if self.mass(under)? >= self.fraction {
        return self.root(below, level);
      }
      let at = self.region(level, LevelBound::At);
      if self.mass(at)? >= self.fraction {
        return self.within_jump(under, at);
      }
      below = level;
    }
    make_internal_error!(
      "Maximum posterior region found no density level that holds the fraction {}",
      self.fraction
    )
  }

  fn region(&self, level: f64, bound: LevelBound) -> Interval {
    Interval {
      lower: self.earlier.position(level, bound),
      upper: self.later.position(level, bound),
    }
  }

  fn mass(&self, region: Interval) -> Result<f64, Report> {
    Ok(self.fractions.mass_between(region.lower, region.upper)? / self.fractions.total())
  }

  fn within_jump(&self, under: Interval, at: Interval) -> Result<Interval, Report> {
    let later_extended = Interval {
      lower: under.lower,
      upper: at.upper,
    };
    if self.mass(later_extended)? >= self.fraction {
      let upper = self.fractions.quantile_from(under.lower, self.fraction);
      return Ok(Interval {
        lower: under.lower,
        upper: upper.clamp(under.upper, at.upper),
      });
    }
    let lower = self.fractions.survival_inverse_from(at.upper, self.fraction);
    Ok(Interval {
      lower: lower.clamp(at.lower, under.lower),
      upper: at.upper,
    })
  }

  fn root(&self, below: f64, above: f64) -> Result<Interval, Report> {
    let above = if above.is_finite() {
      above
    } else {
      self.finite_level_above(below)?
    };
    let problem = LevelMass { solve: self, above };
    let result = Executor::new(problem, BrentRoot::new(below, above, 0.0))
      .configure(|cfg| cfg.max_iters(LEVEL_SOLVE_MAX_ITERS))
      .run()
      .map_err(|e| make_report!("Density level of the maximum posterior region not found: {e}"))?;
    let status = &result.state.termination_status;
    if *status != TerminationStatus::Terminated(TerminationReason::SolverConverged) {
      return make_internal_error!(
        "Density level solve of the maximum posterior region between levels {below} and {above} did not converge: {status:?}"
      );
    }
    let Some(level) = result.state.best_param else {
      return make_internal_error!("Density level solve of the maximum posterior region returned no level");
    };
    Ok(self.region(level, LevelBound::of(level, above)))
  }

  fn finite_level_above(&self, below: f64) -> Result<f64, Report> {
    let mut step = 1.0;
    for _ in 0..LEVEL_BRACKET_DOUBLINGS {
      let above = below + step;
      if self.mass(self.region(above, LevelBound::At))? >= self.fraction {
        return Ok(above);
      }
      step *= 2.0;
    }
    make_internal_error!(
      "Maximum posterior region found no finite density level above {below} that holds the fraction {}",
      self.fraction
    )
  }
}

struct LevelMass<'a> {
  solve: &'a LevelSolve<'a>,
  above: f64,
}

impl CostFunction for LevelMass<'_> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, level: &f64) -> Result<f64, Error> {
    let region = self.solve.region(*level, LevelBound::of(*level, self.above));
    let mass = self
      .solve
      .mass(region)
      .map_err(|report| Error::msg(format!("{report:#}")))?;
    Ok(mass - self.solve.fraction)
  }
}

#[derive(Clone, Copy, Debug)]
struct Interval {
  lower: f64,
  upper: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LevelBound {
  Below,
  At,
}

impl LevelBound {
  fn of(level: f64, above: f64) -> Self {
    if level >= above {
      LevelBound::Below
    } else {
      LevelBound::At
    }
  }
}

struct LevelSide {
  t: Array1<f64>,
  y: Array1<f64>,
  envelope: Array1<f64>,
  side: Side,
  tail: Option<BoundaryBehavior>,
  hard_cell: Option<HardCell>,
  edge_segments: Vec<(usize, EdgeCell)>,
}

impl LevelSide {
  #[expect(
    clippy::needless_pass_by_value,
    reason = "array views are Copy; clippy misses the Copy impl of ndarray 0.17 views"
  )]
  fn new(
    t: ArrayView1<f64>,
    y: ArrayView1<f64>,
    side: Side,
    tail: Option<BoundaryBehavior>,
    hard_cell: Option<HardCell>,
    edge_segments: Vec<(usize, EdgeCell)>,
  ) -> Self {
    let envelope = y
      .iter()
      .scan(f64::NEG_INFINITY, |running_max, &h| {
        *running_max = running_max.max(h);
        Some(*running_max)
      })
      .collect();
    Self {
      t: t.to_owned(),
      y: y.to_owned(),
      envelope,
      side,
      tail,
      hard_cell,
      edge_segments,
    }
  }

  fn position(&self, level: f64, bound: LevelBound) -> f64 {
    let inside = match bound {
      LevelBound::Below => partition_point(&self.envelope, |&h| h < level),
      LevelBound::At => partition_point(&self.envelope, |&h| h <= level),
    };
    if inside == 0 {
      return self.t[0];
    }
    let n_points = self.t.len();
    if inside == n_points {
      return self.tail_position(level);
    }
    if let Some((_, cell)) = self.edge_segments.iter().find(|(segment, _)| *segment == inside) {
      return cell.position_between(level, self.t[inside - 1], self.t[inside]);
    }
    let (y0, y1) = (self.y[inside - 1], self.y[inside]);
    let (x0, x1) = (self.t[inside - 1], self.t[inside]);
    if y1.is_infinite() {
      return self
        .hard_cell
        .map_or(x0, |cell| cell.law.position_at(cell.anchor, level));
    }
    let fraction = (-(level - y0)).exp_m1() / (-(y1 - y0)).exp_m1();
    x0 + fraction.clamp(0.0, 1.0) * (x1 - x0)
  }

  fn tail_position(&self, level: f64) -> f64 {
    let n_points = self.t.len();
    let edge = GridEdge {
      t: self.t[n_points - 1],
      y: self.y[n_points - 1],
    };
    match self.tail {
      Some(BoundaryBehavior::Linear(law)) => law.position_at(edge, level, self.side),
      Some(BoundaryBehavior::HardApproach(law)) => law.position_at(edge, level),
      Some(BoundaryBehavior::Error | BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_)) | None => edge.t,
    }
  }
}
