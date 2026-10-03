use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_primitives::Seq;

#[derive(Debug)]
pub struct OptimizeOutputMaps {
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}

#[derive(Debug, Clone)]
pub struct OptimizeNodeOut {
  pub name: Option<String>,
  pub branch_support: Option<f64>,
}
