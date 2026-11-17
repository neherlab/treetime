use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gamma::{GAMMA_SHAPE_MAX, GAMMA_SHAPE_MIN};
use crate::{make_internal_report, make_report};
use argmin::core::{CostFunction, Error, Executor, State, TerminationReason};
use eyre::Report;
use std::cell::Cell;
use std::f64::consts::GOLDEN_RATIO;
use treetime_primitives::LogLh;

const SHAPE_MAX_ITER: u64 = 100;

const SHAPE_ABS_TOL: f64 = 1e-4;

const SHAPE_STEP: f64 = 0.1;

const SHAPE_STEP_GROWTH: f64 = GOLDEN_RATIO;

pub(crate) fn optimize_gamma_shape(
  start: f64,
  log_lh: impl Fn(f64) -> Result<LogLh, Report>,
) -> Result<ShapeEstimate, Report> {
  let evaluations = Cell::new(0_usize);
  let cost = |log_shape: f64| -> Result<f64, Report> {
    evaluations.set(evaluations.get() + 1);
    Ok(-log_lh(shape_at(log_shape))?.value())
  };
  let log_shape = match bracket_minimum(start.ln(), &cost)? {
    Bracket::Bound(at_bound) => at_bound.x,
    Bracket::Interior { a, b, c } => minimize_in_bracket(a, b, c, &cost)?,
  };
  Ok(ShapeEstimate {
    shape: shape_at(log_shape),
    evaluations: evaluations.get(),
  })
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ShapeEstimate {
  pub(crate) shape: f64,
  pub(crate) evaluations: usize,
}

fn shape_at(log_shape: f64) -> f64 {
  log_shape.exp().clamp(GAMMA_SHAPE_MIN, GAMMA_SHAPE_MAX)
}

fn bracket_minimum(start: f64, cost: &impl Fn(f64) -> Result<f64, Report>) -> Result<Bracket, Report> {
  let start = Point::at(start.clamp(GAMMA_SHAPE_MIN.ln(), GAMMA_SHAPE_MAX.ln()), cost)?;
  let above = neighbor(start, 1.0, cost)?;
  if let Some(above) = above.filter(|above| above.f < start.f) {
    return expand(start, above, cost);
  }
  let below = neighbor(start, -1.0, cost)?;
  if let Some(below) = below.filter(|below| below.f < start.f) {
    return expand(start, below, cost);
  }
  match (below, above) {
    (Some(below), Some(above)) => Ok(Bracket::Interior {
      a: below.x,
      b: start,
      c: above.x,
    }),
    (Some(inner), None) | (None, Some(inner)) => resolve_bound(start, inner.x, cost),
    (None, None) => Ok(Bracket::Bound(start)),
  }
}

#[expect(
  clippy::float_cmp,
  reason = "the clamp returns the start itself when the start lies on the bound in this direction"
)]
fn neighbor(start: Point, direction: f64, cost: &impl Fn(f64) -> Result<f64, Report>) -> Result<Option<Point>, Report> {
  let x = (start.x + direction * SHAPE_STEP).clamp(GAMMA_SHAPE_MIN.ln(), GAMMA_SHAPE_MAX.ln());
  if x == start.x {
    Ok(None)
  } else {
    Ok(Some(Point::at(x, cost)?))
  }
}

#[expect(
  clippy::float_cmp,
  clippy::while_float,
  reason = "the expansion clamps each step to the bound, so it reaches the bound exactly"
)]
fn expand(start: Point, downhill: Point, cost: &impl Fn(f64) -> Result<f64, Report>) -> Result<Bracket, Report> {
  let bound = if downhill.x > start.x {
    GAMMA_SHAPE_MAX.ln()
  } else {
    GAMMA_SHAPE_MIN.ln()
  };
  let (mut previous, mut current) = (start, downhill);
  while current.x != bound {
    let x = current.x + SHAPE_STEP_GROWTH * (current.x - previous.x);
    let next = Point::at(x.clamp(GAMMA_SHAPE_MIN.ln(), GAMMA_SHAPE_MAX.ln()), cost)?;
    if next.f >= current.f {
      return Ok(Bracket::Interior {
        a: previous.x,
        b: current,
        c: next.x,
      });
    }
    (previous, current) = (current, next);
  }
  resolve_bound(current, previous.x, cost)
}

fn resolve_bound(at_bound: Point, inner: f64, cost: &impl Fn(f64) -> Result<f64, Report>) -> Result<Bracket, Report> {
  let offset = SHAPE_ABS_TOL.min((inner - at_bound.x).abs() / 2.0);
  let near = Point::at(at_bound.x + (inner - at_bound.x).signum() * offset, cost)?;
  if at_bound.f <= near.f {
    Ok(Bracket::Bound(at_bound))
  } else {
    Ok(Bracket::Interior {
      a: at_bound.x,
      b: near,
      c: inner,
    })
  }
}

fn minimize_in_bracket(a: f64, b: Point, c: f64, cost: &impl Fn(f64) -> Result<f64, Report>) -> Result<f64, Report> {
  let solver = BrentBracketed::builder()
    .xa(a)
    .xb(b.x)
    .xc(c)
    .fxb(b.f)
    .abs_tol(SHAPE_ABS_TOL)
    .build();
  let shape_cost = ShapeCost { cost };
  let res = Executor::new(&shape_cost, solver)
    .configure(|cfg| cfg.max_iters(SHAPE_MAX_ITER))
    .run()
    .map_err(|e| make_report!("Gamma shape optimization failed: {e}"))?;
  let state = res.state();
  let log_shape = state
    .best_param
    .ok_or_else(|| make_internal_report!("Gamma shape optimization finished without a best parameter"))?;
  match state.get_termination_reason() {
    Some(TerminationReason::SolverConverged | TerminationReason::MaxItersReached) => Ok(log_shape),
    reason => Err(make_internal_report!(
      "Gamma shape optimization stopped for an unexpected reason: {reason:?}"
    )),
  }
}

#[derive(Clone, Copy, Debug)]
struct Point {
  x: f64,
  f: f64,
}

impl Point {
  fn at(x: f64, cost: &impl Fn(f64) -> Result<f64, Report>) -> Result<Self, Report> {
    Ok(Self { x, f: cost(x)? })
  }
}

enum Bracket {
  Bound(Point),
  Interior { a: f64, b: Point, c: f64 },
}

struct ShapeCost<'a, F> {
  cost: &'a F,
}

impl<F> CostFunction for &ShapeCost<'_, F>
where
  F: Fn(f64) -> Result<f64, Report>,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, log_shape: &Self::Param) -> Result<Self::Output, Error> {
    (self.cost)(*log_shape).map_err(|e| Error::msg(format!("{e:#}")))
  }
}
