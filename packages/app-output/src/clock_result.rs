#[derive(Debug, Clone)]
pub struct ClockNodeOut {
  pub name: Option<String>,
  pub div: f64,
  pub time: Option<f64>,
  pub is_outlier: bool,
  pub bad_branch: bool,
}
