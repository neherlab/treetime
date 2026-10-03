use serde::Serialize;
use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_primitives::Seq;

#[derive(Debug, Default)]
pub struct TimetreeOutputMaps {
  pub root_sequence: Option<Seq>,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TimetreeNodeOut {
  pub name: Option<String>,
  pub branch_support: Option<f64>,
  pub time: Option<f64>,
  pub div: f64,
  pub is_outlier: bool,
  pub bad_branch: bool,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct TimetreeEdgeOut {
  pub branch_length: Option<f64>,
  pub date_branch_length: Option<f64>,
}
