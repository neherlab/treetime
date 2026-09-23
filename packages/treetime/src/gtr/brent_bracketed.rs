use argmin::core::{CostFunction, Error, IterState, KV, Problem, Solver, State, TerminationReason};
use bon::bon;

const GOLDEN: f64 = 0.381_966_011_250_105;

pub struct BrentBracketed {
  tol: f64,
  zeps: f64,
  xb: f64,
  a: f64,
  b: f64,
  x: f64,
  w: f64,
  v: f64,
  fx: f64,
  fw: f64,
  fv: f64,
  e: f64,
  d: f64,
}

#[bon]
impl BrentBracketed {
  #[builder]
  pub(crate) fn new(xa: f64, xb: f64, xc: f64) -> Self {
    let (a, b) = if xa < xc { (xa, xc) } else { (xc, xa) };
    Self {
      tol: f64::EPSILON.sqrt(),
      zeps: 1e-11,
      xb,
      a,
      b,
      x: xb,
      w: xb,
      v: xb,
      fx: f64::NAN,
      fw: f64::NAN,
      fv: f64::NAN,
      e: 0.0,
      d: 0.0,
    }
  }
}

impl<O> Solver<O, IterState<f64, (), (), (), (), f64>> for BrentBracketed
where
  O: CostFunction<Param = f64, Output = f64>,
{
  fn name(&self) -> &'static str {
    "BrentBracketed"
  }

  fn init(
    &mut self,
    problem: &mut Problem<O>,
    state: IterState<f64, (), (), (), (), f64>,
  ) -> Result<(IterState<f64, (), (), (), (), f64>, Option<KV>), Error> {
    let fb = problem.cost(&self.xb)?;
    self.fx = fb;
    self.fw = fb;
    self.fv = fb;
    Ok((state.param(self.x).cost(self.fx), None))
  }

  #[expect(
    clippy::many_single_char_names,
    clippy::float_cmp,
    reason = "single-letter names follow the notation of the formulas; Brent's method compares iterates exactly, as the reference algorithm does"
  )]
  fn next_iter(
    &mut self,
    problem: &mut Problem<O>,
    state: IterState<f64, (), (), (), (), f64>,
  ) -> Result<(IterState<f64, (), (), (), (), f64>, Option<KV>), Error> {
    let m = 0.5 * (self.a + self.b);
    let tol1 = self.tol * self.x.abs() + self.zeps;
    let tol2 = 2.0 * tol1;

    if (self.x - m).abs() <= tol2 - 0.5 * (self.b - self.a) {
      return Ok((
        state
          .terminate_with(TerminationReason::SolverConverged)
          .param(self.x)
          .cost(self.fx),
        None,
      ));
    }

    let mut use_golden = true;
    if self.e.abs() > tol1 {
      let r = (self.x - self.w) * (self.fx - self.fv);
      let q0 = (self.x - self.v) * (self.fx - self.fw);
      let mut p = (self.x - self.v) * q0 - (self.x - self.w) * r;
      let mut q = 2.0 * (q0 - r);
      if q > 0.0 {
        p = -p;
      } else {
        q = -q;
      }
      let old_e = self.e;
      self.e = self.d;

      if p.abs() < (0.5 * q * old_e).abs() && p > q * (self.a - self.x) && p < q * (self.b - self.x) {
        self.d = p / q;
        let u = self.x + self.d;
        if (u - self.a) < tol2 || (self.b - u) < tol2 {
          self.d = if m > self.x { tol1 } else { -tol1 };
        }
        use_golden = false;
      }
    }

    if use_golden {
      self.e = if self.x < m { self.b - self.x } else { self.a - self.x };
      self.d = GOLDEN * self.e;
    }

    let u = if self.d.abs() >= tol1 {
      self.x + self.d
    } else {
      self.x + if self.d > 0.0 { tol1 } else { -tol1 }
    };
    let fu = problem.cost(&u)?;

    if fu <= self.fx {
      if u < self.x {
        self.b = self.x;
      } else {
        self.a = self.x;
      }
      self.v = self.w;
      self.fv = self.fw;
      self.w = self.x;
      self.fw = self.fx;
      self.x = u;
      self.fx = fu;
    } else {
      if u < self.x {
        self.a = u;
      } else {
        self.b = u;
      }
      if fu <= self.fw || self.w == self.x {
        self.v = self.w;
        self.fv = self.fw;
        self.w = u;
        self.fw = fu;
      } else if fu <= self.fv || self.v == self.x || self.v == self.w {
        self.v = u;
        self.fv = fu;
      }
    }

    Ok((state.param(self.x).cost(self.fx), None))
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use argmin::core::Executor;
  use eyre::Report;

  struct Fn1D(fn(f64) -> f64);

  impl CostFunction for Fn1D {
    type Param = f64;
    type Output = f64;

    fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
      Ok((self.0)(*x))
    }
  }

  fn minimize(f: fn(f64) -> f64, xa: f64, xb: f64, xc: f64) -> Result<f64, Report> {
    let solver = BrentBracketed::builder().xa(xa).xb(xb).xc(xc).build();
    let res = Executor::new(Fn1D(f), solver)
      .configure(|cfg| cfg.max_iters(500))
      .run()
      .map_err(|e| eyre::eyre!("{e}"))?;
    res
      .state()
      .best_param
      .ok_or_else(|| eyre::eyre!("solver reported no best_param"))
  }

  #[test]
  fn test_brent_bracketed_quadratic() -> Result<(), Report> {
    let xmin = minimize(|x| (x - 3.0) * (x - 3.0), 0.0, 3.5, 10.0)?;
    assert_abs_diff_eq!(xmin, 3.0, epsilon = 1e-7);
    Ok(())
  }

  #[test]
  fn test_brent_bracketed_shifted_quadratic() -> Result<(), Report> {
    let xmin = minimize(|x| (x - 0.7) * (x - 0.7) + 1.0, 0.0, 1.0, 5.0)?;
    assert_abs_diff_eq!(xmin, 0.7, epsilon = 1e-7);
    Ok(())
  }

  #[test]
  fn test_brent_bracketed_cosine() -> Result<(), Report> {
    let xmin = minimize(f64::cos, 2.0, 3.0, 5.0)?;
    assert_abs_diff_eq!(xmin, std::f64::consts::PI, epsilon = 1e-7);
    Ok(())
  }

  #[test]
  fn test_brent_bracketed_seeds_at_interior_point() -> Result<(), Report> {
    let xmin = minimize(|x| (x - 5.0) * (x - 5.0), 0.0, 4.9, 10.0)?;
    assert_abs_diff_eq!(xmin, 5.0, epsilon = 1e-7);
    Ok(())
  }
}
