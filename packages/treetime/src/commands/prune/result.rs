use crate::gtr::gtr::GTR;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::mutation::Mutation;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[derive(Serialize)]
pub struct PruneGraphData {
  pub gtr: Option<GTR>,
}

impl PruneGraphData {
  pub fn new(gtr: Option<GTR>) -> Self {
    Self { gtr }
  }
}

/// Nucleotide sequences and mutations gathered from the prune partition for the tree writers.
///
/// Gathered once, serially, from the partition while it is in scope in the command, so the auspice,
/// phyloxml, and MAT writers read plain value maps instead of reading the partition during
/// serialization.
#[derive(Debug, Default)]
pub struct PruneOutputMaps {
  /// Reconstructed nucleotide root sequence, or `None` when no partition exists.
  pub root_sequence: Option<Seq>,
  /// Reconstructed nucleotide sequence per node, for phyloxml clades.
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
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
/// name/support and per-edge branch length the writers consume, and `gtr` holds the fitted model.
/// `graph` carries the tree the output writers still read from.
#[derive(Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: GraphAncestral<PruneGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, PruneNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
}
