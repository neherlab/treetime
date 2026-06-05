use ndarray::Array1;

/// Shared storage and validation for piecewise functions on non-uniform breakpoints.
///
/// Holds breakpoints (strictly ascending) and associated values. The relationship
/// between the number of breakpoints and values depends on the function type:
/// piecewise constant uses `n+1` values (one per region), piecewise linear uses
/// `n` values (one per breakpoint).
#[derive(Debug, Clone)]
pub struct PiecewiseFnBase {
  breakpoints: Array1<f64>,
  values: Array1<f64>,
}

impl PiecewiseFnBase {
  pub fn new(breakpoints: Array1<f64>, values: Array1<f64>) -> Self {
    debug_assert!(
      breakpoints
        .as_slice()
        .unwrap()
        .windows(2)
        .all(|w| matches!(w, [a, b] if a < b)),
      "breakpoints must be strictly ascending"
    );
    Self { breakpoints, values }
  }

  pub fn breakpoints(&self) -> &Array1<f64> {
    &self.breakpoints
  }

  pub fn values(&self) -> &Array1<f64> {
    &self.values
  }

  pub fn breakpoints_slice(&self) -> &[f64] {
    self.breakpoints.as_slice().unwrap()
  }

  pub fn values_slice(&self) -> &[f64] {
    self.values.as_slice().unwrap()
  }
}
