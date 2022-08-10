use argmin::core::observers::{ObserverMode, SlogLogger};
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;
use eyre::{eyre, Report};
use log::log_enabled;
use log::Level::Trace;

struct CostFunctionWrapper<F>
where
  F: Fn(f64) -> f64,
{
  problem: F,
}

impl<F> CostFunctionWrapper<F>
where
  F: Fn(f64) -> f64,
{
  pub const fn new(problem: F) -> Self {
    Self { problem }
  }
}

impl<F> CostFunction for CostFunctionWrapper<F>
where
  F: Fn(f64) -> f64,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
    let problem = &self.problem;
    Ok(problem(*x))
  }
}

pub fn minimize_scalar_brent_bounded(problem: impl Fn(f64) -> f64, bounds: (f64, f64)) -> Result<(f64, f64), Report> {
  let problem = CostFunctionWrapper::new(problem);
  let solver = BrentOpt::new(bounds.0, bounds.1);

  let mut executor = Executor::new(problem, solver).configure(|state| state.max_iters(1000));

  if log_enabled!(Trace) {
    executor = executor.add_observer(SlogLogger::term_noblock(), ObserverMode::NewBest);
  }

  let result = executor.run().map_err(|err| eyre!("{err}"))?;

  let param = *result
    .state()
    .get_best_param()
    .ok_or_else(|| eyre!("Unable to get the best param"))?;

  let cost = result.state().get_best_cost();

  Ok((param, cost))
}
