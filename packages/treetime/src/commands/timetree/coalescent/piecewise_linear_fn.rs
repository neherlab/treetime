#![allow(dead_code)]
use ndarray::Array1;

/// Piecewise linear function represented by breakpoints and values.
///
/// For n breakpoints, there are n values (one per breakpoint).
/// Evaluation uses linear interpolation between breakpoints.
/// Extrapolation uses constant values at boundaries.
///
/// Breakpoints must be sorted in ascending order.
#[derive(Debug, Clone)]
pub struct PiecewiseLinearFn {
  breakpoints: Array1<f64>,
  values: Array1<f64>,
}

impl PiecewiseLinearFn {
  /// Create from breakpoints and values.
  ///
  /// # Arguments
  /// - `breakpoints`: sorted ascending, n elements (n >= 2)
  /// - `values`: n elements, one per breakpoint
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
    debug_assert!(
      breakpoints
        .as_slice()
        .unwrap()
        .windows(2)
        .all(|w| matches!(w, [a, b] if a < b))
    );
    Self { breakpoints, values }
  }

  /// Get breakpoint times.
  pub fn breakpoints(&self) -> &Array1<f64> {
    &self.breakpoints
  }

  /// Get values at breakpoints.
  pub fn values(&self) -> &Array1<f64> {
    &self.values
  }

  /// Evaluate at a single point using linear interpolation.
  pub fn eval(&self, t: f64) -> f64 {
    let n = self.breakpoints.len();

    if t <= self.breakpoints[0] {
      return self.values[0];
    }
    if t >= self.breakpoints[n - 1] {
      return self.values[n - 1];
    }

    let idx = self.breakpoints.as_slice().unwrap().partition_point(|&bp| bp < t);
    let i = idx.saturating_sub(1);

    let t0 = self.breakpoints[i];
    let t1 = self.breakpoints[i + 1];
    let y0 = self.values[i];
    let y1 = self.values[i + 1];

    let alpha = (t - t0) / (t1 - t0);
    y0 + alpha * (y1 - y0)
  }

  /// Evaluate at multiple points (must be monotonically sorted ascending).
  pub fn eval_many(&self, queries: &Array1<f64>) -> Array1<f64> {
    debug_assert!(queries.windows(2).into_iter().all(|w| w[0] <= w[1]));

    let n = self.breakpoints.len();
    let breakpoints = self.breakpoints.as_slice().unwrap();
    let values = self.values.as_slice().unwrap();
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_piecewise_linear_eval_interpolation() {
    // Breakpoints: [0.0, 10.0]
    // Values: [0.0, 100.0]
    // Linear from 0 to 100
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);

    pretty_assert_abs_diff_eq!(pl.eval(0.0), 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(5.0), 50.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(10.0), 100.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(2.5), 25.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(7.5), 75.0, epsilon = 1e-9);
  }

  #[test]
  fn test_piecewise_linear_eval_extrapolation() {
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);

    // Before first breakpoint: constant at first value
    pretty_assert_abs_diff_eq!(pl.eval(-5.0), 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(-100.0), 0.0, epsilon = 1e-9);

    // After last breakpoint: constant at last value
    pretty_assert_abs_diff_eq!(pl.eval(15.0), 100.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(1000.0), 100.0, epsilon = 1e-9);
  }

  #[test]
  fn test_piecewise_linear_eval_multiple_segments() {
    // Three segments: 0->10: 0->50, 10->20: 50->100, 20->30: 100->0
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0, 20.0, 30.0], array![0.0, 50.0, 100.0, 0.0]);

    // First segment
    pretty_assert_abs_diff_eq!(pl.eval(5.0), 25.0, epsilon = 1e-9);

    // Second segment
    pretty_assert_abs_diff_eq!(pl.eval(15.0), 75.0, epsilon = 1e-9);

    // Third segment (decreasing)
    pretty_assert_abs_diff_eq!(pl.eval(25.0), 50.0, epsilon = 1e-9);

    // Exact breakpoints
    pretty_assert_abs_diff_eq!(pl.eval(10.0), 50.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pl.eval(20.0), 100.0, epsilon = 1e-9);
  }

  #[test]
  fn test_piecewise_linear_eval_many() {
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);
    let ts = array![-5.0, 0.0, 5.0, 10.0, 15.0];
    let result = pl.eval_many(&ts);

    pretty_assert_abs_diff_eq!(result[0], 0.0, epsilon = 1e-9); // extrapolation before
    pretty_assert_abs_diff_eq!(result[1], 0.0, epsilon = 1e-9); // at first breakpoint
    pretty_assert_abs_diff_eq!(result[2], 50.0, epsilon = 1e-9); // interpolation
    pretty_assert_abs_diff_eq!(result[3], 100.0, epsilon = 1e-9); // at last breakpoint
    pretty_assert_abs_diff_eq!(result[4], 100.0, epsilon = 1e-9); // extrapolation after
  }

  #[test]
  fn test_piecewise_linear_accessors() {
    let breakpoints = array![1.0, 5.0, 10.0];
    let values = array![10.0, 20.0, 30.0];
    let pl = PiecewiseLinearFn::new(breakpoints.clone(), values.clone());

    assert_eq!(pl.breakpoints(), &breakpoints);
    assert_eq!(pl.values(), &values);
  }
}
