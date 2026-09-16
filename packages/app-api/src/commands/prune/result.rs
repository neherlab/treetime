use treetime::seq::mutation::Mutation;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

/// Nucleotide root sequence and mutations gathered from the prune partition for the tree writers.
///
/// Gathered once, serially, from the partition while it is in scope in the command, so the auspice
/// and MAT writers read plain value maps instead of reading the partition during
/// serialization.
#[derive(Debug, Default)]
pub struct PruneOutputMaps {
  /// Reconstructed nucleotide root sequence, or `None` when no partition exists.
  pub root_sequence: Option<Seq>,
  /// Nucleotide mutations (substitutions followed by indels) per edge.
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

/// Per-node prune output as a value: the name and input branch support the output writers read.
#[derive(Debug, Clone, Serialize)]
pub struct PruneNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

/// Per-edge prune output as a value: the branch length the output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}

/// Prune result as a value.
///
/// The durable outputs are reachable directly off the result: `nodes` and `edges` hold the per-node
/// name/support and per-edge branch length the writers consume. `graph` carries the tree topology the
/// output writers read from.
#[derive(Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, PruneNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}
