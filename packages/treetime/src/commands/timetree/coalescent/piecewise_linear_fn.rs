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
  #[cfg_attr(not(test), allow(dead_code))]
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
  #[cfg_attr(not(test), allow(dead_code))]
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
