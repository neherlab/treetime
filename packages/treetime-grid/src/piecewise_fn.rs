use ndarray::Array1;

#[derive(Debug, Clone)]
pub struct PiecewiseFnBase {
  breakpoints: Array1<f64>,
  values: Array1<f64>,
}

impl PiecewiseFnBase {
  pub(crate) fn new(breakpoints: Array1<f64>, values: Array1<f64>) -> Self {
    debug_assert!(
      breakpoints.iter().is_sorted_by(|a, b| a < b),
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
}
