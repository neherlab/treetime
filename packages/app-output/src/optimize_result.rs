use serde::Serialize;
use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[derive(Debug)]
pub struct OptimizeOutputMaps {
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}

pub struct OptimizeResult {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, OptimizeNodeOut>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

#[derive(Debug, Clone, Serialize)]
pub struct OptimizeNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}
