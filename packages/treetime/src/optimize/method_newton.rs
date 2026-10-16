use crate::optimize::indel::IndelTerm;
use crate::optimize::likelihood::evaluate_with_indels;
use crate::optimize::solver_outcome::{Convergence, EdgeSolution};
use crate::optimize::zero_boundary::search_upper_bound;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use num::clamp;

const NEWTON_REL_TOL: f64 = 0.001;

const NEWTON_ABS_TOL: f64 = 1e-8;

const NEWTON_MAX_ITER: usize = 50;

pub(super) fn newton_sqrt_inner(
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
  let (interval, mut current) = initial_point(start, interval, &evaluate)?;
  let mut best = current;
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
      return Ok(best.converged());
    }
    if best_is_current && bracket.width() <= tolerance {
      return Ok(best.converged());
    }
    current = Evaluated::at(proposal.s(), &evaluate)?;
    best_is_current = current.point.log_lh >= best.point.log_lh;
    if best_is_current {
      best = current;
    }
  }

  Ok(SqrtSolution {
    s: best.s,
    convergence: Convergence::IterationLimit,
  })
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

  fn including(self, s: f64) -> Self {
    Self {
      lower: self.lower.min(s),
      upper: self.upper.max(s),
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

fn initial_point(
  start: f64,
  interval: SqrtInterval,
  evaluate: impl Fn(f64) -> Result<SqrtPoint, Report>,
) -> Result<(SqrtInterval, Evaluated), Report> {
  if interval.contains(start) {
    return Ok((interval, Evaluated::at(start, evaluate)?));
  }
  let at_start = Evaluated::at(start, &evaluate)?;
  if start >= 0.0 && at_start.point.log_lh.is_finite() {
    Ok((interval.including(start), at_start))
  } else {
    let clamped = clamp(start, interval.lower, interval.upper);
    Ok((interval, Evaluated::at(clamped, evaluate)?))
  }
}

fn is_concave_maximum(point: &SqrtPoint, tolerance: f64) -> bool {
  point.curvature <= 0.0 && point.slope.abs() <= tolerance * point.curvature.abs()
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

  fn converged(self) -> SqrtSolution {
    SqrtSolution {
      s: self.s,
      convergence: Convergence::Converged,
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
  f64::max(0.5 * NEWTON_REL_TOL * s, NEWTON_ABS_TOL)
}

pub(super) fn sqrt_step_lower_bound(s: f64) -> f64 {
  s - (s * s + 1.0).sqrt()
}

pub(super) fn chain_rule_sqrt(s: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_ds = 2.0 * s * dl_dt;
  let d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
  (dl_ds, d2l_ds2)
}
