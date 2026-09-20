use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct ClockNodeOut {
  pub name: Option<String>,
  pub div: f64,
  pub time: Option<f64>,
  pub is_outlier: bool,
  pub bad_branch: bool,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}
