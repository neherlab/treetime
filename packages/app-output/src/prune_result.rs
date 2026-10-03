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

pub struct PruneResult {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, PruneNodeOut>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

#[derive(Debug, Clone)]
pub struct PruneNodeOut {
  pub name: Option<String>,
  pub branch_support: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}
