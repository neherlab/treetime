use serde::Serialize;
use std::collections::BTreeMap;
use treetime::seq::mutation::{Mutation, Sub};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Debug, Default)]
pub struct AncestralOutputMaps {
  pub root_sequence: Option<Seq>,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

#[derive(Debug)]
pub struct AugurOutputMaps {
  pub root_sequence: Seq,
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  pub edge_subs: BTreeMap<GraphEdgeKey, Vec<Sub>>,
  pub sequence_length: usize,
  pub ambiguous_char: AsciiChar,
}

#[derive(Serialize)]
pub struct AncestralResult {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, AncestralNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

#[derive(Debug, Clone, Serialize)]
pub struct AncestralNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}
