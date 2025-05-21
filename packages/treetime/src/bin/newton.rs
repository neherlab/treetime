/// Calculate square root of a number by solving equation `x^2 - n = 0` using Newton-Raphson method.
///
/// The Newton solver in `argmin` requires `Gradient` and `Hessian` traits to be implemented for the problem.
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{Error, Executor, Gradient, Hessian, KV, State};
use argmin::solver::newton::Newton;
use ctor::ctor;
use log::{LevelFilter, info};
use std::fmt::{Debug, Display};
use treetime::utils::global_init::{global_init, setup_logger};

fn main() -> Result<(), Box<dyn std::error::Error>> {
  let n = 2.0;
  let problem = SqrtSearch { n };
  let solver = Newton::new().with_gamma(1.0)?;
  let res = Executor::new(problem, solver)
    .configure(move |cfg| cfg.param(n).max_iters(100).target_cost(1e-8))
    // .add_observer(SlogLogger::term_noblock(), ObserverMode::Always)
    .add_observer(MyObserver {}, ObserverMode::Always)
    .run()?
    .state
    .best_param
    .unwrap();

  println!("The square root of {n} is approximately {res}");
  Ok(())
}

/// Problem definition
struct SqrtSearch {
  // Current state of the parameter
  n: f64,
}

// /// Cost function: `(x^2 - n)^2` - we are aiming to minimize this to find the square root.
// /// `CostFunction` trait is not required for `Newton` solver, so it is commented away here,
// /// but can be used with other solvers
// impl CostFunction for SqrtFunction {
//   type Param = f64;
//   type Output = f64;
//
//   fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
//     Ok((x.powi(2) - self.n).powi(2))
//   }
// }

/// Gradient of the cost function: derivative with respect to x, which is `4x(x^2 - n)`.
/// This represents the first derivative used in Newton's method to find the direction of update.
impl Gradient for SqrtSearch {
  type Param = f64;
  type Gradient = f64;

  fn gradient(&self, x: &Self::Param) -> Result<Self::Gradient, Error> {
    Ok(4.0 * x * (x.powi(2) - self.n))
  }
}

// Hessian of the cost function: second derivative with respect to x, which is `12x^2 - 4n`.
// This represents the second derivative used in Newton's method to adjust the step size optimally.
impl Hessian for SqrtSearch {
  type Param = f64;
  type Hessian = f64;

  fn hessian(&self, x: &Self::Param) -> Result<Self::Hessian, Error> {
    Ok(12.0 * x.powi(2) - 4.0 * self.n)
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
      "param={:.9}, best_param={:.9}",
      state.get_param().unwrap(),
      state.get_best_param().unwrap()
    );
    Ok(())
  }
}

#[ctor]
fn init() {
  global_init();
  setup_logger(LevelFilter::Info);
}
