use serde::Serialize;
use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[derive(Debug, Default)]
pub struct PruneOutputMaps {
  pub root_sequence: Option<Seq>,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

#[derive(Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, PruneNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PruneNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}
