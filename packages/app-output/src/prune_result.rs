use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_primitives::Seq;

#[derive(Debug, Default)]
pub struct PruneOutputMaps {
  pub root_sequence: Option<Seq>,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

#[derive(Debug, Clone)]
pub struct PruneNodeOut {
  pub name: Option<String>,
  pub branch_support: Option<f64>,
}
