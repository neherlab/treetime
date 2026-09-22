use ndarray::Array1;

#[derive(Debug, Clone)]
pub struct PiecewiseFnBase {
  breakpoints: Array1<f64>,
  values: Array1<f64>,
}

impl PiecewiseFnBase {
  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  pub(crate) fn new(breakpoints: Array1<f64>, values: Array1<f64>) -> Self {
    debug_assert!(
      breakpoints.as_slice().unwrap().is_sorted_by(|a, b| a < b),
      "breakpoints must be strictly ascending"
    );
    Self { breakpoints, values }
  }

  pub(crate) fn breakpoints(&self) -> &Array1<f64> {
    &self.breakpoints
  }

  pub(crate) fn values(&self) -> &Array1<f64> {
    &self.values
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  pub(crate) fn breakpoints_slice(&self) -> &[f64] {
    self.breakpoints.as_slice().unwrap()
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  pub(crate) fn values_slice(&self) -> &[f64] {
    self.values.as_slice().unwrap()
  }
}
