//! One-dimensional minimization by Brent's method

use crate::{make_internal_report, make_report};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;
use eyre::Report;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BrentMinimum {
  pub x: f64,
  pub value: f64,
}

/// Minimizes `f` on `[lower, upper]` in at most `max_iters` iterations. `tolerance` is argmin's `(relative, absolute)`
/// pair; `None` keeps argmin's default.
pub(crate) fn minimize_brent(
  lower: f64,
  upper: f64,
  tolerance: Option<(f64, f64)>,
  max_iters: u64,
  f: impl Fn(f64) -> Result<f64, Report>,
) -> Result<BrentMinimum, Report> {
  let solver = BrentOpt::new(lower, upper);
  let solver = match tolerance {
    Some((relative, absolute)) => solver.set_tolerance(relative, absolute),
    None => solver,
  };
  let result = Executor::new(ClosureCost(f), solver)
    .configure(|cfg| cfg.max_iters(max_iters))
    .run()
    .map_err(|error| make_report!("Brent minimization on [{lower}, {upper}] failed: {error}"))?;
  let state = result.state();
  let x = state
    .best_param
    .ok_or_else(|| make_internal_report!("Brent minimization on [{lower}, {upper}] reported no minimum"))?;
  Ok(BrentMinimum {
    x,
    value: state.best_cost,
  })
}

struct ClosureCost<F>(F);

impl<F: Fn(f64) -> Result<f64, Report>> CostFunction for ClosureCost<F> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
    (self.0)(*x).map_err(|error| Error::msg(error.to_string()))
  }
}
