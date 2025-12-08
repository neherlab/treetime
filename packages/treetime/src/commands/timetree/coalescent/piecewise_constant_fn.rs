#![allow(dead_code)]
use ndarray::Array1;

/// Piecewise constant function represented by breakpoints and values.
///
/// For n breakpoints, there are n+1 value regions:
/// - values[0] for t < breakpoints[0]
/// - values[i] for breakpoints[i-1] <= t < breakpoints[i]  (1 <= i < n)
/// - values[n] for t >= breakpoints[n-1]
///
/// Breakpoints must be sorted in ascending order.
#[derive(Debug, Clone)]
pub struct PiecewiseConstantFn {
  breakpoints: Array1<f64>,
  values: Array1<f64>,
}

impl PiecewiseConstantFn {
  /// Create from breakpoints and values.
  ///
  /// # Arguments
  /// - `breakpoints`: sorted ascending, n elements
  /// - `values`: n+1 elements
  pub fn new(breakpoints: Array1<f64>, values: Array1<f64>) -> Self {
    debug_assert!(breakpoints.len() + 1 == values.len());
    debug_assert!(
      breakpoints
        .as_slice()
        .unwrap()
        .windows(2)
        .all(|w| matches!(w, [a, b] if a < b))
    );
    Self { breakpoints, values }
  }

  /// Get breakpoint times (where function changes value).
  pub fn breakpoints(&self) -> &Array1<f64> {
    &self.breakpoints
  }

  /// Evaluate at a single point.
  pub fn eval(&self, t: f64) -> f64 {
    let idx = self.breakpoints.as_slice().unwrap().partition_point(|&bp| bp <= t);
    self.values[idx]
  }

  /// Evaluate at multiple points (must be monotonically sorted ascending).
  pub fn eval_many(&self, queries: &Array1<f64>) -> Array1<f64> {
    debug_assert!(queries.windows(2).into_iter().all(|w| w[0] <= w[1]));

    let breakpoints = self.breakpoints.as_slice().unwrap();
    let mut bp_iter = breakpoints.iter().peekable();

    queries.mapv(|t| {
      while bp_iter.peek().is_some_and(|&&bp| bp <= t) {
        bp_iter.next();
      }
      self.values[breakpoints.len() - bp_iter.len()]
    })
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_piecewise_constant_eval() {
    // Breakpoints: [1.0, 5.0, 10.0]
    // Values: [0.0, 1.0, 2.0, 3.0]
    // t < 1.0 -> 0.0
    // 1.0 <= t < 5.0 -> 1.0
    // 5.0 <= t < 10.0 -> 2.0
    // t >= 10.0 -> 3.0
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_abs_diff_eq!(pc.eval(-1.0), 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(0.5), 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(1.0), 1.0, epsilon = 1e-9); // at breakpoint: value after
    pretty_assert_abs_diff_eq!(pc.eval(3.0), 1.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(5.0), 2.0, epsilon = 1e-9); // at breakpoint: value after
    pretty_assert_abs_diff_eq!(pc.eval(7.5), 2.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(10.0), 3.0, epsilon = 1e-9); // at breakpoint: value after
    pretty_assert_abs_diff_eq!(pc.eval(100.0), 3.0, epsilon = 1e-9);
  }

  #[test]
  fn test_piecewise_constant_eval_many() {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0], array![0.0, 1.0, 2.0]);
    let ts = array![0.0, 1.0, 3.0, 5.0, 10.0];
    let result = pc.eval_many(&ts);
    pretty_assert_abs_diff_eq!(result[0], 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(result[1], 1.0, epsilon = 1e-9); // at breakpoint: value after
    pretty_assert_abs_diff_eq!(result[2], 1.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(result[3], 2.0, epsilon = 1e-9); // at breakpoint: value after
    pretty_assert_abs_diff_eq!(result[4], 2.0, epsilon = 1e-9);
  }
}
