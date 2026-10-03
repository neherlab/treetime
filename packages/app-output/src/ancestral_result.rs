use std::collections::BTreeMap;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Debug)]
pub struct AncestralOutputMaps {
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

#[derive(Debug)]
pub struct AugurOutputMaps {
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  pub sequence_length: usize,
  pub ambiguous_char: AsciiChar,
}

#[derive(Debug, Clone)]
pub struct AncestralNodeOut {
  pub name: Option<String>,
  pub branch_support: Option<f64>,
}
