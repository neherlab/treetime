use serde::Serialize;

pub trait TestCase: Clone + Send + Sync + Serialize {
  fn name(&self) -> &str;

  fn description(&self) -> &str;

  fn stress_type(&self) -> &str;

  fn analytical_caution(&self) -> &str;

  fn dx(&self) -> f64;
}
