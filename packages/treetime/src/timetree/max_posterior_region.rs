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
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
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
  trim: f64,
) -> Result<Option<PosteriorSurvival>, Report> {
  match Posterior::of(distribution)? {
    Posterior::Empty => Ok(None),
    Posterior::Point(t) => Ok(Some(PosteriorSurvival {
      t: vec![t],
      survival: vec![0.0],
    })),
    Posterior::Range(start, end) => Ok(Some(PosteriorSurvival {
      t: vec![start, end],
      survival: vec![1.0, 0.0],
    })),
    Posterior::Function(function) => function_survival(&function, trim).map(Some),
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

fn function_survival(function: &DistributionFunction<f64, NegLog>, trim: f64) -> Result<PosteriorSurvival, Report> {
  let profile = MassProfile::new(function)?;
  let fractions = posterior_fractions(&profile)?;
  let (x_min, x_max) = (profile.x_min(), profile.x_max());
  let (dx, n_points) = (function.dx(), function.y().len());
  let earlier = tail_points(x_min - fractions.quantile(trim), dx, n_points).map(|distance| x_min - distance);
  let grid = (0..n_points).map(|idx| x_min + dx * idx.to_f64_exact());
  let later = tail_points(fractions.survival_inverse(trim) - x_max, dx, n_points).map(|distance| x_max + distance);
  let t = chain!(earlier.rev(), grid, later).collect_vec();
  let survival = t.iter().map(|&x| fractions.survival(x)).collect();
  Ok(PosteriorSurvival { t, survival })
}

fn tail_points(distance: f64, dx: f64, max_points: usize) -> impl DoubleEndedIterator<Item = f64> {
  let n_points = if distance.is_finite() && distance > 0.0 {
    f64_to_usize_saturating((distance / dx).ceil()).clamp(1, max_points)
  } else {
    0
  };
  (1..=n_points).map(move |k| distance * (k.to_f64_exact() / n_points.to_f64_exact()))
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
