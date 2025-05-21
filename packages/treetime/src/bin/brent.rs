/// Calculate π by finding a root of the cosine function.
/// Brent's method is used to find the root of `cos(x) = 0` which corresponds to `x = π/2`.
/// This implementation seeks to find π by solving for π/2 and doubling the result.
///
/// The Brent solver in `argmin` requires `CostFunction` trait to be implemented for the problem.
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{CostFunction, Error, Executor, KV, State};
use argmin::solver::brent::BrentOpt;
use ctor::ctor;
use log::{LevelFilter, info};
use std::fmt::{Debug, Display};
use treetime::utils::global_init::{global_init, setup_logger};

fn main() -> Result<(), Box<dyn std::error::Error>> {
  let problem = PiSearch;
  let bounds = 1.0..2.0; // π/2 is ~1.57
  let solver = BrentOpt::new(bounds.start, bounds.end);

  let res = Executor::new(problem, solver)
    .configure(move |cfg| cfg.max_iters(100).target_cost(1e-8))
    // .add_observer(SlogLogger::term_noblock(), ObserverMode::Always)
    .add_observer(MyObserver {}, ObserverMode::Always)
    .run()?
    .state
    .best_param
    .unwrap();

  println!("The value of π is approximately {}", res * 2.0);
  Ok(())
}

/// Problem definition
struct PiSearch;

/// Cost function for finding Pi.
/// The cosine function crosses zero at `π/2`, `3π/2`, `5π/2`, etc.
/// Here, we aim to find the first positive root, i.e., `π/2`.
/// This function returns the square of the cosine of x, `cos(x)^2`, which has a minimum at `π/2`.
/// We use the square to ensure a non-negative output that reaches zero at our root of interest.
impl CostFunction for PiSearch {
  type Param = f64;
  type Output = f64;

  fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
    Ok(x.cos().powi(2))
  }
}

/// Custom observer which prints optimization state on every iteration
#[derive(Debug)]
struct MyObserver;

impl<I> Observe<I> for MyObserver
where
  I: State<Param: Display> + Debug,
{
  // Is executed after initialization of a solver
  fn observe_init(&mut self, name: &str, state: &I, kv: &KV) -> Result<(), Error> {
    info!("{name}");
    self.observe_iter(state, kv)
  }

  // Is executed after each iteration of a solver
  fn observe_iter(&mut self, state: &I, _kv: &KV) -> Result<(), Error> {
    info!(
      "iter={:>06}, param={:.9}, best_param={:.9}, cost={:.9}, best_cost={:.9}",
      state.get_iter(),
      state.get_param().unwrap(),
      state.get_best_param().unwrap(),
      state.get_cost(),
      state.get_best_cost(),
    );
    Ok(())
  }
}

#[ctor]
fn init() {
  global_init();
  setup_logger(LevelFilter::Info);
}
