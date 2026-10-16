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
  let mut s = clamp(start, interval.lower, interval.upper);
  let mut point = evaluate(s)?;
  let mut best = BestPoint {
    s,
    log_lh: point.log_lh,
  };
  let mut bracket = interval;

  for _ in 0..max_iter {
    bracket = bracket.narrowed(s, point.slope);
    let next = newton_step(s, &point, &interval, &bracket).unwrap_or_else(|| bracket.midpoint());
    let converged = (next - s).abs() <= newton_tolerance_sqrt(s);
    point = evaluate(next)?;
    if point.log_lh > best.log_lh {
      best = BestPoint {
        s: next,
        log_lh: point.log_lh,
      };
    }
    if converged {
      return Ok(SqrtSolution {
        s: best.s,
        convergence: Convergence::Converged,
      });
    }
    s = next;
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

  fn contains(&self, s: f64) -> bool {
    (self.lower..=self.upper).contains(&s)
  }

  fn midpoint(&self) -> f64 {
    f64::midpoint(self.lower, self.upper)
  }
}

fn newton_step(s: f64, point: &SqrtPoint, interval: &SqrtInterval, bracket: &SqrtInterval) -> Option<f64> {
  let curvature = -point.curvature.abs();
  let stationary_minimum = point.slope == 0.0 && point.curvature > 0.0;
  (curvature < 0.0 && !stationary_minimum)
    .then(|| {
      let step = clamp(point.slope / curvature, sqrt_step_lower_bound(s), s);
      clamp(s - step, interval.lower, interval.upper)
    })
    .filter(|&next| bracket.contains(next))
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

struct BestPoint {
  s: f64,
  log_lh: f64,
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
