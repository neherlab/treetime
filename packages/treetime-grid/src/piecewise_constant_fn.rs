use crate::piecewise_fn::PiecewiseFnBase;
use itertools::Itertools;
use ndarray::Array1;

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
    let idx = self.base.breakpoints_slice().partition_point(|&bp| bp <= t);
    self.base.values()[idx]
  }

  #[cfg_attr(
    dylint_lib = "treetime_lints",
    allow(
      pub_unused_in_workspace,
      reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach"
    )
  )]
  pub fn eval_left(&self, t: f64) -> f64 {
    let idx = self.base.breakpoints_slice().partition_point(|&bp| bp < t);
    self.base.values()[idx]
  }

  #[must_use]
  pub fn zip_map(&self, other: &Self, f: impl Fn(f64, f64) -> f64) -> Self {
    let mut breakpoints = self
      .breakpoints()
      .iter()
      .chain(other.breakpoints())
      .copied()
      .sorted_by(f64::total_cmp)
      .collect_vec();
    breakpoints.dedup_by(|left, right| left.total_cmp(right).is_eq());

    let values = std::iter::once(f(self.values()[0], other.values()[0]))
      .chain(breakpoints.iter().map(|&t| f(self.eval(t), other.eval(t))))
      .collect();

    Self::new(Array1::from(breakpoints), values)
  }
}
