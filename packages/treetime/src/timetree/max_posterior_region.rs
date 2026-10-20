use argmin::core::{CostFunction, Error as ArgminError, Executor, TerminationReason, TerminationStatus};
use argmin::solver::brent::BrentRoot;
use eyre::Report;
use itertools::{Itertools, chain};
use ndarray::{Array1, ArrayView1, s};
use ordered_float::OrderedFloat;
use std::iter::once;
use treetime_distribution::{
  BoundaryBehavior, Distribution, DistributionFunction, MassFractions, MassProfile, NegLog, TailMass,
};
use treetime_grid::{GridEdge, HardApproachLaw, HardCell, Side};
use treetime_utils::array::ndarray::partition_point;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error, make_report};

const LEVEL_BRACKET_DOUBLINGS: usize = 64;
const LEVEL_SOLVE_MAX_ITERS: u64 = 200;
const REGION_MASS_ROUNDING_ULPS: f64 = 64.0;
const HERMITE_SLOPE_BASIS_MAX: f64 = 4.0 / 27.0;

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
  trim: f64,
  max_error: f64,
) -> Result<Option<PosteriorSurvival>, Report> {
  match Posterior::of(distribution)? {
    Posterior::Empty => Ok(None),
    Posterior::Point(t) => Ok(Some(PosteriorSurvival::point(t))),
    Posterior::Range(start, end) => {
      let slope = -(end - start).recip();
      Ok(Some(PosteriorSurvival {
        t: vec![start, end],
        survival: vec![1.0, 0.0],
        slope: vec![slope, slope],
      }))
    },
    Posterior::Function(function) => function_survival(&function, trim, max_error).map(Some),
  }
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
  pub slope: Vec<f64>,
}

impl PosteriorSurvival {
  fn point(t: f64) -> Self {
    Self {
      t: vec![t],
      survival: vec![0.0],
      slope: vec![0.0],
    }
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
  trim: f64,
  max_error: f64,
) -> Result<PosteriorSurvival, Report> {
  let profile = MassProfile::new(function)?;
  let fractions = posterior_fractions(&profile)?;
  let start = fractions.quantile(trim);
  let end = fractions.survival_inverse(trim).max(start);
  if !(start.is_finite() && end.is_finite()) {
    let (side, tail) = if start.is_finite() {
      (Side::Right, function.right_extrap())
    } else {
      (Side::Left, function.left_extrap())
    };
    return match tail {
      BoundaryBehavior::Linear(law) => make_error!(
        "Posterior on [{}, {}] has a tail of slope {:e} on the {side:?} side, which holds mass {trim:e} beyond the largest finite time",
        function.x_min(),
        function.x_max(),
        law.slope
      ),
      BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) | BoundaryBehavior::Error => {
        make_internal_error!("Posterior quantile {trim} on the {side:?} side without a soft tail is not finite")
      },
    };
  }
  let pieces = survival_pieces(function, start, end)?;
  if pieces.is_empty() {
    return Ok(PosteriorSurvival::point(start));
  }
  let mut knots = SurvivalKnots {
    fractions: &fractions,
    max_error,
    survival: PosteriorSurvival {
      t: vec![],
      survival: vec![],
      slope: vec![],
    },
  };
  for piece in pieces {
    knots.add_piece(piece);
  }
  Ok(knots.survival)
}

fn survival_pieces(
  function: &DistributionFunction<f64, NegLog>,
  start: f64,
  end: f64,
) -> Result<Vec<SurvivalPiece>, Report> {
  let y = function.y();
  let n = y.len();
  let x_at = |idx: usize| function.x_min() + function.dx() * idx.to_f64_exact();
  let (x_min, x_max) = (x_at(0), x_at(n - 1));
  let left = SurvivalPiece {
    start,
    end: end.min(x_min),
    law: PieceLaw::of_tail(function.left_extrap(), Side::Left),
  };
  let right = SurvivalPiece {
    start: start.max(x_max),
    end,
    law: PieceLaw::of_tail(function.right_extrap(), Side::Right),
  };
  let cells = (0..n - 1)
    .filter(|&cell| x_at(cell + 1) > start && x_at(cell) < end)
    .map(|cell| {
      let law = match (y[cell].is_infinite(), y[cell + 1].is_infinite()) {
        (false, false) | (true, true) => PieceLaw::Exact,
        (hard_left, _) => PieceLaw::PowerLaw {
          t_hard: x_at(if hard_left { cell } else { cell + 1 }),
          b: function.hard_cell(cell)?.map_or(0.0, |hard| hard.law.b),
        },
      };
      Ok(SurvivalPiece {
        start: start.max(x_at(cell)),
        end: end.min(x_at(cell + 1)),
        law,
      })
    });
  chain!(once(Ok(left)), cells, once(Ok(right)))
    .filter_ok(|piece| piece.start < piece.end)
    .try_collect()
}

#[derive(Clone, Copy)]
struct SurvivalPiece {
  start: f64,
  end: f64,
  law: PieceLaw,
}

#[derive(Clone, Copy)]
enum PieceLaw {
  Exact,
  Exponential { rate: f64 },
  PowerLaw { t_hard: f64, b: f64 },
}

impl PieceLaw {
  fn of_tail(behavior: BoundaryBehavior, side: Side) -> Self {
    match (behavior, side) {
      (BoundaryBehavior::Linear(law), Side::Left) => Self::Exponential { rate: -law.slope },
      (BoundaryBehavior::Linear(law), Side::Right) => Self::Exponential { rate: law.slope },
      (BoundaryBehavior::HardApproach(law), _) => Self::PowerLaw {
        t_hard: law.t_hard,
        b: law.b,
      },
      (BoundaryBehavior::Hard | BoundaryBehavior::Error, _) => Self::Exact,
    }
  }

  fn interpolation_error_bound(self, a: Knot, b: Knot) -> f64 {
    let h = b.t - a.t;
    let quartic = h.powi(4) / 384.0;
    let smooth = match self {
      PieceLaw::Exact => 0.0,
      PieceLaw::Exponential { rate } => rate.abs().powi(3) * a.slope.abs().max(b.slope.abs()) * quartic,
      PieceLaw::PowerLaw { t_hard, b: exponent } => {
        let factor = (exponent * (exponent - 1.0) * (exponent - 2.0)).abs();
        let (near, far) = if (a.t - t_hard).abs() <= (b.t - t_hard).abs() {
          (a, b)
        } else {
          (b, a)
        };
        let largest_at = if exponent < 3.0 { near } else { far };
        let distance = (largest_at.t - t_hard).abs();
        if factor <= 0.0 {
          0.0
        } else if distance > 0.0 {
          factor * largest_at.slope.abs() / distance.powi(3) * quartic
        } else {
          f64::INFINITY
        }
      },
    };
    let monotone = (a.survival - b.survival).abs() + HERMITE_SLOPE_BASIS_MAX * h * (a.slope.abs() + b.slope.abs());
    if smooth.is_finite() {
      smooth.min(monotone)
    } else {
      monotone
    }
  }
}

#[derive(Clone, Copy)]
struct Knot {
  t: f64,
  survival: f64,
  slope: f64,
}

struct SurvivalKnots<'a> {
  fractions: &'a MassFractions<'a>,
  max_error: f64,
  survival: PosteriorSurvival,
}

impl SurvivalKnots<'_> {
  fn add_piece(&mut self, piece: SurvivalPiece) {
    let a = self.knot(piece.start, Side::Right);
    let b = self.knot(piece.end, Side::Left);
    self.push(a);
    self.refine(piece.law, a, b);
  }

  fn knot(&self, t: f64, side: Side) -> Knot {
    Knot {
      t,
      survival: self.fractions.survival(t),
      slope: -self.fractions.density(t, side),
    }
  }

  fn refine(&mut self, law: PieceLaw, a: Knot, b: Knot) {
    let t_mid = f64::midpoint(a.t, b.t);
    if !(a.t < t_mid && t_mid < b.t) || law.interpolation_error_bound(a, b) <= self.max_error {
      self.push(b);
      return;
    }
    let mid = self.knot(t_mid, Side::Left);
    self.refine(law, a, mid);
    self.refine(law, mid, b);
  }

  fn push(&mut self, knot: Knot) {
    let PosteriorSurvival { t, survival, slope } = &mut self.survival;
    let repeated = t.last().is_some_and(|&last| last.total_cmp(&knot.t).is_eq())
      && slope.last().is_some_and(|&last| last.total_cmp(&knot.slope).is_eq());
    if !repeated {
      t.push(knot.t);
      survival.push(knot.survival);
      slope.push(knot.slope);
    }
  }
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
  let y = function.y();
  let n_points = y.len();
  let x_at = |idx: usize| function.x_min() + function.dx() * idx.to_f64_exact();
  let (x_min, x_max) = (x_at(0), x_at(n_points - 1));
  let left_carries_mass = fractions.cdf(x_min) > 0.0;
  let right_carries_mass = fractions.survival(x_max) > 0.0;
  if peak == n_points - 1 && !right_carries_mass {
    return Ok(PosteriorRegion {
      lower: fractions.survival_inverse(fraction),
      upper: x_max,
      peak: PeakLocation::Latest,
    });
  }
  if peak == 0 && !left_carries_mass {
    return Ok(PosteriorRegion {
      lower: x_min,
      upper: fractions.quantile(fraction),
      peak: PeakLocation::Earliest,
    });
  }

  let t = Array1::from_shape_fn(n_points, x_at);
  let earlier = LevelSide::new(
    t.slice(s![..=peak;-1]),
    y.slice(s![..=peak;-1]),
    Side::Left,
    left_carries_mass.then_some(function.left_extrap()),
    first_hard_cell(function, (0..peak).rev(), |cell| cell)?,
  );
  let later = LevelSide::new(
    t.slice(s![peak..]),
    y.slice(s![peak..]),
    Side::Right,
    right_carries_mass.then_some(function.right_extrap()),
    first_hard_cell(function, peak..n_points - 1, |cell| cell + 1)?,
  );
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

  fn cost(&self, level: &f64) -> Result<f64, ArgminError> {
    let region = self.solve.region(*level, LevelBound::of(*level, self.above));
    let mass = self
      .solve
      .mass(region)
      .map_err(|report| ArgminError::msg(format!("{report:#}")))?;
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
      Some(BoundaryBehavior::Error | BoundaryBehavior::Hard) | None => edge.t,
    }
  }
}
