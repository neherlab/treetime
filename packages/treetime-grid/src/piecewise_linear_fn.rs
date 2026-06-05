use crate::piecewise_fn::PiecewiseFnBase;
use ndarray::Array1;

/// Piecewise linear function on non-uniform breakpoints with constant extrapolation.
///
/// Represents a continuous function f: R -> R defined by n breakpoints and n values
/// (one value per breakpoint). Between consecutive breakpoints the function
/// interpolates linearly. Outside the breakpoint range the function extends at the
/// boundary values (constant extrapolation):
///
/// ```text
///   f(t) = values[0]                                    for t <= breakpoints[0]
///   f(t) = lerp(values[i], values[i+1], alpha)          for breakpoints[i] <= t <= breakpoints[i+1]
///   f(t) = values[n-1]                                  for t >= breakpoints[n-1]
/// ```
///
/// where `alpha = (t - breakpoints[i]) / (breakpoints[i+1] - breakpoints[i])`.
///
/// This is the non-uniform counterpart of [`GridFn`](crate::GridFn), which
/// represents the same kind of function on a uniform grid and uses index
/// arithmetic for O(1) point location. `PiecewiseLinearFn` imposes no spacing
/// constraint and uses binary search for point location in O(log n).
///
/// # Invariants
///
/// - `breakpoints` is strictly ascending (enforced by debug assertion)
/// - `breakpoints.len() == values.len()`
/// - At least 2 breakpoints
#[derive(Debug, Clone)]
pub struct PiecewiseLinearFn {
  base: PiecewiseFnBase,
}

impl PiecewiseLinearFn {
  /// Create from breakpoints and values.
  ///
  /// # Arguments
  /// - `breakpoints`: sorted strictly ascending, n elements (n >= 2)
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
    Self {
      base: PiecewiseFnBase::new(breakpoints, values),
    }
  }

  /// Breakpoint positions where the slope may change.
  pub fn breakpoints(&self) -> &Array1<f64> {
    self.base.breakpoints()
  }

  /// Values at breakpoints.
  #[allow(dead_code)]
  pub fn values(&self) -> &Array1<f64> {
    self.base.values()
  }

  /// Evaluate at a single point using linear interpolation.
  ///
  /// Returns the linearly interpolated value within the enclosing segment,
  /// or the nearest boundary value when `t` falls outside the breakpoint range.
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

  /// Evaluate at multiple points in a single sweep.
  ///
  /// Takes advantage of sorted queries to walk the breakpoint array in tandem,
  /// achieving O(n + m) total time instead of O(m log n) for independent lookups.
  ///
  /// # Precondition
  ///
  /// `queries` must be sorted in non-decreasing order (debug-asserted).
  #[allow(dead_code)]
  pub fn eval_many(&self, queries: &Array1<f64>) -> Array1<f64> {
    debug_assert!(queries.windows(2).into_iter().all(|w| w[0] <= w[1]));

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
