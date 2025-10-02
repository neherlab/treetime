use crate::commands::timetree::branch_distributions::BranchDistribution;
use crate::commands::timetree::date_constraints::DateConstraint;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::node::GraphNodeKey;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct NodeTemporalData {
  pub constraint: Option<DateConstraint>,
  pub bad_branch: bool,
  pub time_before_present: Option<f64>,
  pub posterior: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct EdgeTemporalData {
  pub clock_length: Option<f64>,
  pub distribution: Option<BranchDistribution>,
  pub msg_to_parent: Option<Arc<Distribution>>,
  pub msg_to_child: Option<Arc<Distribution>>,
  pub msg_from_child: Option<Arc<Distribution>>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TemporalDataStore {
  pub nodes: BTreeMap<GraphNodeKey, NodeTemporalData>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeTemporalData>,
  pub sequence_length: Option<usize>,
  pub iteration: Option<usize>,
}
