use serde::{Deserialize, Serialize};

pub trait TestCase: Clone + Send + Sync + Serialize {
  fn base(&self) -> &TestCaseBase;

  fn name(&self) -> &str {
    self.base().name()
  }

  fn description(&self) -> &str {
    self.base().description()
  }

  fn stress_type(&self) -> &str {
    self.base().stress_type()
  }

  fn analytical_caution(&self) -> &str {
    self.base().analytical_caution()
  }

  fn slowness(&self) -> f64 {
    self.base().slowness()
  }

  fn input_grid_domain(&self) -> (f64, f64);

  fn input_grid_n_points(&self) -> usize;
}

/// Metadata fields shared by every test case across all test suites.
///
/// Embedded in each concrete test case with `#[serde(flatten)]` so the
/// serialized form stays flat (fields are not nested under a `base` key). The
/// accessors here are the single implementation forwarded to by the `TestCase`
/// trait defaults.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestCaseBase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
}

impl TestCaseBase {
  pub fn name(&self) -> &str {
    &self.name
  }

  pub fn description(&self) -> &str {
    &self.description
  }

  pub fn stress_type(&self) -> &str {
    &self.stress_type
  }

  pub fn analytical_caution(&self) -> &str {
    &self.analytical_caution
  }

  pub fn slowness(&self) -> f64 {
    self.slowness
  }
}
