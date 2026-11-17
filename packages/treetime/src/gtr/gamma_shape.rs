use crate::gtr::gamma::{GAMMA_SHAPE_MAX, GAMMA_SHAPE_MIN};
use crate::{make_internal_report, make_report};
use argmin::core::{CostFunction, Error, Executor, State, TerminationReason};
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use treetime_primitives::LogLh;

const SHAPE_MAX_ITER: u64 = 100;

const SHAPE_REL_TOL: f64 = 1.490_116_119_384_765_6e-8;

const SHAPE_ABS_TOL: f64 = 1e-6;

pub(crate) fn optimize_gamma_shape(log_lh: impl Fn(f64) -> Result<LogLh, Report>) -> Result<f64, Report> {
  let cost_fn = ShapeCost { log_lh: &log_lh };
  let solver = BrentOpt::new(GAMMA_SHAPE_MIN.ln(), GAMMA_SHAPE_MAX.ln()).set_tolerance(SHAPE_REL_TOL, SHAPE_ABS_TOL);
  let res = Executor::new(&cost_fn, solver)
    .configure(|cfg| cfg.max_iters(SHAPE_MAX_ITER))
    .run()
    .map_err(|e| make_report!("Gamma shape optimization failed: {e}"))?;
  let state = res.state();
  let log_shape = state
    .best_param
    .ok_or_else(|| make_internal_report!("Gamma shape optimization finished without a best parameter"))?;
  match state.get_termination_reason() {
    Some(TerminationReason::SolverConverged | TerminationReason::MaxItersReached) => Ok(shape_at(log_shape)),
    reason => Err(make_internal_report!(
      "Gamma shape optimization stopped for an unexpected reason: {reason:?}"
    )),
  }
}

fn shape_at(log_shape: f64) -> f64 {
  log_shape.exp().clamp(GAMMA_SHAPE_MIN, GAMMA_SHAPE_MAX)
}

struct ShapeCost<'a, F> {
  log_lh: &'a F,
}

impl<F> CostFunction for &ShapeCost<'_, F>
where
  F: Fn(f64) -> Result<LogLh, Report>,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, log_shape: &Self::Param) -> Result<Self::Output, Error> {
    let log_lh = (self.log_lh)(shape_at(*log_shape)).map_err(|e| Error::msg(format!("{e:#}")))?;
    Ok(-log_lh.value())
  }
}
