use ndarray::Array1;

/// Piecewise constant (step) function on non-uniform breakpoints.
///
/// Represents a right-continuous step function f: R -> R defined by n breakpoints
/// and n+1 values. The function is constant between consecutive breakpoints and
/// jumps to a new value at each breakpoint:
///
/// ```text
///   f(t) = values[0]           for t < breakpoints[0]
///   f(t) = values[i]           for breakpoints[i-1] <= t < breakpoints[i]   (1 <= i < n)
///   f(t) = values[n]           for t >= breakpoints[n-1]
/// ```
///
/// Breakpoints must be strictly ascending. The function is right-continuous at each
/// breakpoint: `eval(bp)` returns the value in the interval starting at `bp`.
/// For left-continuous evaluation (the value just before the breakpoint), use
/// [`eval_left`](Self::eval_left).
///
/// This is the non-uniform counterpart of [`GridFn`](crate::GridFn), which
/// represents a piecewise linear function on a uniform grid. `PiecewiseConstantFn`
/// requires no uniform spacing constraint and uses binary search for point location
/// instead of index arithmetic.
///
/// # Invariants
///
/// - `breakpoints` is strictly ascending (enforced by debug assertion)
/// - `values.len() == breakpoints.len() + 1`
#[derive(Debug, Clone)]
pub struct PiecewiseConstantFn {
  breakpoints: Array1<f64>,
  values: Array1<f64>,
}

impl PiecewiseConstantFn {
  /// Create from breakpoints and values.
  ///
  /// # Arguments
  /// - `breakpoints`: sorted strictly ascending, n elements
  /// - `values`: n+1 elements (one per inter-breakpoint region, plus the two tails)
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

  /// Breakpoint positions where the function changes value.
  pub fn breakpoints(&self) -> &Array1<f64> {
    &self.breakpoints
  }

  /// Evaluate at a single point (right-continuous).
  ///
  /// At a breakpoint, returns the value of the interval that starts at that
  /// breakpoint (the post-jump value).
  pub fn eval(&self, t: f64) -> f64 {
    let idx = self.breakpoints.as_slice().unwrap().partition_point(|&bp| bp <= t);
    self.values[idx]
  }

  /// Evaluate the left limit at a single point.
  ///
  /// At a breakpoint, returns the value of the interval ending at that
  /// breakpoint (the pre-jump value). Between breakpoints and away from
  /// breakpoints, behaves identically to [`eval`](Self::eval).
  pub fn eval_left(&self, t: f64) -> f64 {
    let idx = self.breakpoints.as_slice().unwrap().partition_point(|&bp| bp < t);
    self.values[idx]
  }

  /// Evaluate at multiple points in a single sweep.
  ///
  /// Takes advantage of sorted queries to avoid repeated binary searches:
  /// walks the breakpoint array in tandem with the query array in O(n + m)
  /// time, where n is the number of breakpoints and m is the number of queries.
  ///
  /// # Precondition
  ///
  /// `queries` must be sorted in non-decreasing order (debug-asserted).
  #[allow(dead_code)]
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
