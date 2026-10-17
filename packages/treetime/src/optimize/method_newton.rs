use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels;
use crate::optimize::solver_outcome::{Convergence, EdgeSolution};
use crate::optimize::zero_boundary::search_upper_bound;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use itertools::Itertools;
use num::clamp;
use ordered_float::OrderedFloat;

const NEWTON_REL_TOL: f64 = 0.001;

const NEWTON_ABS_TOL: f64 = 1e-8;

const NEWTON_MAX_BRANCH_LENGTH_ERROR: f64 = 1e-6;

const NEWTON_MAX_ITER: usize = 50;

pub(super) fn solve_newton(
  branch_length: f64,
  contributions: &[OptimizationContribution],
  indels: IndelTerm,
  min_branch_length: f64,
  one_mutation: f64,
) -> Result<EdgeSolution, Report> {
  let interval = SqrtInterval {
    lower: min_branch_length.sqrt(),
    upper: search_upper_bound(branch_length, one_mutation).sqrt(),
  };
  let solution = newton_sqrt_maximize(branch_length.sqrt(), interval, NEWTON_MAX_ITER, |s| {
    let metrics = evaluate_with_indels(contributions, indels, s * s)?;
    let (slope, curvature) = chain_rule_sqrt(s, metrics.derivative, metrics.second_derivative);
    Ok(SqrtPoint {
      log_lh: metrics.log_lh.value(),
      slope,
      curvature,
    })
  })?;
  Ok(EdgeSolution {
    branch_length: solution.s * solution.s,
    convergence: solution.convergence,
  })
}

pub(super) fn newton_sqrt_maximize(
  start: f64,
  interval: SqrtInterval,
  max_iter: usize,
  evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>,
) -> Result<SqrtSolution, Report> {
  let incumbent = feasible_start_outside(start, &interval, &evaluate)?;
  let first = Evaluated::at(clamp(start, interval.lower, interval.upper), &evaluate)?;
  let search = search(first, interval, max_iter, &evaluate)?;
  let best = match incumbent {
    Some(incumbent) if incumbent.point.log_lh > search.best.point.log_lh => incumbent,
    _ => search.best,
  };
  Ok(SqrtSolution {
    s: best.s,
    convergence: search.convergence,
  })
}

fn feasible_start_outside(
  start: f64,
  interval: &SqrtInterval,
  evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>,
) -> Result<Option<Evaluated>, Report> {
  if interval.contains(start) {
    return Ok(None);
  }
  let at_start = Evaluated::at(start, evaluate)?;
  Ok((start >= 0.0 && at_start.point.log_lh.is_finite()).then_some(at_start))
}

fn search(
  first: Evaluated,
  interval: SqrtInterval,
  max_iter: usize,
  evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>,
) -> Result<Search, Report> {
  let mut current = first;
  let mut best = first;
  let mut best_is_current = true;
  let mut bracket = interval;

  for _ in 0..max_iter {
    bracket = bracket.narrowed(current.s, current.point.slope);
    let tolerance = newton_tolerance_sqrt(current.s);
    let proposal = propose(&current, &interval, &bracket, tolerance);
    if best_is_current && is_concave_maximum(&current.point, tolerance) {
      if let Proposal::Newton(next) = proposal {
        best = best.better(Evaluated::at(next, &evaluate)?);
      }
      return finish(best, &interval, &evaluate);
    }
    if best_is_current && is_degenerate_stationary(&current.point) {
      match higher_neighbor(&current, &interval, tolerance, &evaluate)? {
        None => return finish(best, &interval, &evaluate),
        Some(neighbor) => {
          best = neighbor;
          best_is_current = false;
        },
      }
    }
    if best_is_current && bracket.width() <= tolerance {
      return finish(best, &interval, &evaluate);
    }
    current = Evaluated::at(proposal.s(), &evaluate)?;
    best_is_current = current.point.log_lh >= best.point.log_lh;
    if best_is_current {
      best = current;
    }
  }

  Ok(Search {
    best,
    convergence: Convergence::IterationLimit,
  })
}

fn finish(
  best: Evaluated,
  interval: &SqrtInterval,
  evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>,
) -> Result<Search, Report> {
  if interval.upper - best.s > newton_tolerance_sqrt(best.s) {
    return Ok(Search {
      best,
      convergence: Convergence::Converged,
    });
  }
  let at_upper = Evaluated::at(interval.upper, evaluate)?;
  let convergence = if at_upper.point.slope > 0.0 {
    Convergence::SearchBound
  } else {
    Convergence::Converged
  };
  let best = if at_upper.point.log_lh > best.point.log_lh {
    at_upper
  } else {
    best
  };
  Ok(Search { best, convergence })
}

struct Search {
  best: Evaluated,
  convergence: Convergence,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct SqrtInterval {
  pub lower: f64,
  pub upper: f64,
}

impl SqrtInterval {
  fn narrowed(self, s: f64, slope: f64) -> Self {
    if slope > 0.0 {
      Self { lower: s, ..self }
    } else if slope < 0.0 {
      Self { upper: s, ..self }
    } else {
      self
    }
  }

  fn contains(&self, s: f64) -> bool {
    (self.lower..=self.upper).contains(&s)
  }

  fn width(&self) -> f64 {
    self.upper - self.lower
  }

  fn bisection_away_from(&self, s: f64) -> f64 {
    if s - self.lower >= self.upper - s {
      f64::midpoint(self.lower, s)
    } else {
      f64::midpoint(s, self.upper)
    }
  }
}

fn is_concave_maximum(point: &SqrtPoint, tolerance: f64) -> bool {
  point.curvature < 0.0 && point.slope.abs() <= tolerance * point.curvature.abs()
}

fn is_degenerate_stationary(point: &SqrtPoint) -> bool {
  point.slope == 0.0 && point.curvature == 0.0
}

fn higher_neighbor(
  current: &Evaluated,
  interval: &SqrtInterval,
  tolerance: f64,
  evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>,
) -> Result<Option<Evaluated>, Report> {
  let below = (current.s > interval.lower).then(|| f64::max(current.s - tolerance, interval.lower));
  let above = (current.s < interval.upper).then(|| f64::min(current.s + tolerance, interval.upper));
  let neighbors: Vec<Evaluated> = [below, above]
    .into_iter()
    .flatten()
    .map(|s| Evaluated::at(s, &evaluate))
    .try_collect()?;
  Ok(
    neighbors
      .into_iter()
      .filter(|neighbor| neighbor.point.log_lh > current.point.log_lh)
      .max_by_key(|neighbor| OrderedFloat(neighbor.point.log_lh)),
  )
}

fn propose(current: &Evaluated, interval: &SqrtInterval, bracket: &SqrtInterval, tolerance: f64) -> Proposal {
  newton_step(current.s, &current.point, interval, bracket, tolerance).map_or_else(
    || Proposal::Bisection(bracket.bisection_away_from(current.s)),
    Proposal::Newton,
  )
}

fn newton_step(
  s: f64,
  point: &SqrtPoint,
  interval: &SqrtInterval,
  bracket: &SqrtInterval,
  tolerance: f64,
) -> Option<f64> {
  let curvature = -point.curvature.abs();
  let stationary_minimum = point.slope == 0.0 && point.curvature > 0.0;
  (curvature < 0.0 && !stationary_minimum)
    .then(|| {
      let step = clamp(point.slope / curvature, sqrt_step_lower_bound(s), s);
      clamp(s - step, interval.lower, interval.upper)
    })
    .filter(|&next| bracket.contains(next))
    .filter(|&next| point.curvature < 0.0 || (next - s).abs() > tolerance)
}

#[derive(Clone, Copy, Debug)]
enum Proposal {
  Newton(f64),
  Bisection(f64),
}

impl Proposal {
  fn s(self) -> f64 {
    match self {
      Self::Newton(s) | Self::Bisection(s) => s,
    }
  }
}

#[derive(Clone, Copy, Debug)]
struct Evaluated {
  s: f64,
  point: SqrtPoint,
}

impl Evaluated {
  fn at(s: f64, evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>) -> Result<Self, Report> {
    Ok(Self { s, point: evaluate(s)? })
  }

  fn better(self, other: Self) -> Self {
    if other.point.log_lh >= self.point.log_lh {
      other
    } else {
      self
    }
  }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct SqrtPoint {
  pub log_lh: f64,
  pub slope: f64,
  pub curvature: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct SqrtSolution {
  pub s: f64,
  pub convergence: Convergence,
}

pub(super) fn newton_tolerance_sqrt(s: f64) -> f64 {
  let max_branch_length_error = NEWTON_MAX_BRANCH_LENGTH_ERROR / ((s * s + NEWTON_MAX_BRANCH_LENGTH_ERROR).sqrt() + s);
  f64::max(
    f64::min(0.5 * NEWTON_REL_TOL * s, max_branch_length_error),
    NEWTON_ABS_TOL,
  )
}

pub(super) fn sqrt_step_lower_bound(s: f64) -> f64 {
  s - (s * s + 1.0).sqrt()
}

pub(super) fn chain_rule_sqrt(s: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_ds = 2.0 * s * dl_dt;
  let d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
  (dl_ds, d2l_ds2)
}
