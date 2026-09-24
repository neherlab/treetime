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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestCaseBase {
  pub(crate) name: String,
  pub(crate) description: String,
  pub(crate) stress_type: String,
  pub(crate) analytical_caution: String,
  pub(crate) slowness: f64,
}

impl TestCaseBase {
  fn name(&self) -> &str {
    &self.name
  }

  fn description(&self) -> &str {
    &self.description
  }

  fn stress_type(&self) -> &str {
    &self.stress_type
  }

  fn analytical_caution(&self) -> &str {
    &self.analytical_caution
  }

  fn slowness(&self) -> f64 {
    self.slowness
  }
}
