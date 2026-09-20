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

  pub fn breakpoints(&self) -> &Array1<f64> {
    self.base.breakpoints()
  }

  #[allow(dead_code)]
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

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  #[allow(dead_code)]
  pub fn eval_many(&self, queries: &Array1<f64>) -> Array1<f64> {
    debug_assert!(queries.as_slice().unwrap().is_sorted());

    let n = self.base.breakpoints().len();
    let breakpoints = self.base.breakpoints_slice();
    let values = self.base.values_slice();
    let mut bp_iter = breakpoints.iter().enumerate().peekable();

    queries.mapv(|t| {
      if t <= breakpoints[0] {
        return values[0];
      }
      if t >= breakpoints[n - 1] {
        return values[n - 1];
      }

      while bp_iter.peek().is_some_and(|&(_, &bp)| bp < t) {
        bp_iter.next();
      }

      let seg_idx = bp_iter.peek().map_or(n - 2, |(i, _)| i.saturating_sub(1));
      let t0 = breakpoints[seg_idx];
      let t1 = breakpoints[seg_idx + 1];
      let y0 = values[seg_idx];
      let y1 = values[seg_idx + 1];

      let alpha = (t - t0) / (t1 - t0);
      y0 + alpha * (y1 - y0)
    })
  }
}
