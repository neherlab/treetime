use serde::Serialize;

/// Per-node clock output as a value.
///
/// Holds the durable per-node results the clock output writers consume: the estimated `time`
/// (numerical date), the cumulative divergence `div`, and the two exclusion flags. `name` is carried
/// for writers that key by name. The values are gathered from the tree after clock estimation and
/// rerooting complete, so the map is keyed by the final (post-reroot) node set.
#[derive(Debug, Clone, Serialize)]
pub struct ClockNodeOut {
  pub name: Option<String>,
  pub div: f64,
  pub time: Option<f64>,
  pub is_outlier: bool,
  pub bad_branch: bool,
}

/// Per-edge output as a value: the branch length the output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}
