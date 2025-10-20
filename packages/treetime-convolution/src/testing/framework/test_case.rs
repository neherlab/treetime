use serde::Serialize;

pub trait TestCase: Clone + Send + Sync + Serialize {
  fn name(&self) -> &str;

  fn description(&self) -> &str;

  fn stress_type(&self) -> &str;

  fn analytical_caution(&self) -> &str;

  fn input_grid_domain(&self) -> (f64, f64);

  fn input_grid_n_points(&self) -> usize;

  fn output_grid_domain(&self) -> (f64, f64);

  fn output_grid_n_points(&self) -> usize;
}
