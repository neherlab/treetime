use crate::piecewise_fn::PiecewiseFnBase;
use ndarray::Array1;
use treetime_utils::array::ndarray::partition_point;

#[derive(Debug, Clone)]
pub struct PiecewiseConstantFn {
  base: PiecewiseFnBase,
}

impl PiecewiseConstantFn {
  pub fn new(breakpoints: Array1<f64>, values: Array1<f64>) -> Self {
    debug_assert_eq!(breakpoints.len() + 1, values.len());
    Self {
      base: PiecewiseFnBase::new(breakpoints, values),
    }
  }

  pub fn breakpoints(&self) -> &Array1<f64> {
    self.base.breakpoints()
  }

  pub fn values(&self) -> &Array1<f64> {
    self.base.values()
  }

  pub fn eval(&self, t: f64) -> f64 {
    let idx = partition_point(self.base.breakpoints(), |&bp| bp <= t);
    self.base.values()[idx]
  }

  #[cfg_attr(
    dylint_lib = "custom",
    allow(
      pub_unused_in_workspace,
      reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach"
    )
  )]
  pub fn eval_left(&self, t: f64) -> f64 {
    let idx = partition_point(self.base.breakpoints(), |&bp| bp < t);
    self.base.values()[idx]
  }
}
