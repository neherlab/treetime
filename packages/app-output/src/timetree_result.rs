use serde::Serialize;
use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[derive(Debug, Default)]
pub struct TimetreeOutputMaps {
  pub root_sequence: Option<Seq>,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

#[derive(Serialize)]
pub struct TimetreeResult {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TimetreeNodeOut {
  pub name: Option<String>,
  pub desc: Option<String>,
  pub confidence: Option<f64>,
  pub time: Option<f64>,
  pub div: f64,
  pub is_outlier: bool,
  pub bad_branch: bool,
  pub rate_susceptibility_dates: Option<[f64; 3]>,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct TimetreeEdgeOut {
  pub branch_length: Option<f64>,
  pub time_length: Option<f64>,
  pub clock_branch_length: Option<f64>,
  pub gamma: f64,
}

impl TimetreeEdgeOut {
  pub fn profile_branch_length(&self) -> Option<f64> {
    self.clock_branch_length.or(self.branch_length)
  }
}
