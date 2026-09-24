use crate::piecewise_fn::PiecewiseFnBase;
use ndarray::Array1;

#[derive(Debug, Clone)]
pub struct PiecewiseLinearFn {
  base: PiecewiseFnBase,
}

impl PiecewiseLinearFn {
  pub fn new(breakpoints: Array1<f64>, values: Array1<f64>) -> Self {
    debug_assert!(
      breakpoints.len() >= 2,
      "PiecewiseLinearFn requires at least 2 breakpoints"
    );
    debug_assert_eq!(
      breakpoints.len(),
      values.len(),
      "breakpoints and values must have equal length"
    );
    Self {
      base: PiecewiseFnBase::new(breakpoints, values),
    }
  }

  pub fn values(&self) -> &Array1<f64> {
    self.base.values()
  }

  pub fn eval(&self, t: f64) -> f64 {
    let n = self.base.breakpoints().len();

    if t <= self.base.breakpoints()[0] {
      return self.base.values()[0];
    }
    if t >= self.base.breakpoints()[n - 1] {
      return self.base.values()[n - 1];
    }

    let idx = self.base.breakpoints_slice().partition_point(|&bp| bp < t);
    let i = idx.saturating_sub(1);

    let t0 = self.base.breakpoints()[i];
    let t1 = self.base.breakpoints()[i + 1];
    let y0 = self.base.values()[i];
    let y1 = self.base.values()[i + 1];

    let alpha = (t - t0) / (t1 - t0);
    y0 + alpha * (y1 - y0)
  }
}
