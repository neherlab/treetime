use crate::commands::ancestral::aa_node_data::AaNodeData;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::seq::mutation::{Mutation, Sub};
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Serialize)]
pub struct AncestralGraphData {
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  pub mask: Vec<bool>,
  pub aa_node_data: Option<AaNodeData>,
}

impl AncestralGraphData {
  pub fn new(gtr: Option<GTR>, model_name: GtrModelName, mask: Vec<bool>, aa_node_data: Option<AaNodeData>) -> Self {
    Self {
      gtr,
      model_name,
      mask,
      aa_node_data,
    }
  }
}

/// Nucleotide sequences and mutations gathered from the ancestral partition for the tree writers.
///
/// Gathered once, serially, from the partition while it is in scope in the command, so the auspice,
/// phyloxml, MAT, and Newick-comment writers read plain value maps instead of reading the partition.
/// The amino-acid tracks are still merged by the writers from the graph data slot; this struct carries
/// only the nucleotide reads that moved off the partition.
#[derive(Debug, Default)]
pub struct AncestralOutputMaps {
  /// Reconstructed nucleotide root sequence (posterior-resolved), or `None` when no partition exists.
  pub root_sequence: Option<Seq>,
  /// Reconstructed nucleotide sequence per node (posterior-resolved argmax), for phyloxml clades.
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  /// Nucleotide mutations (substitutions followed by indels) per edge.
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

/// Sequences and substitutions gathered from the ancestral partition for the augur node-data writer.
///
/// The augur writer resolves the dense per-node sequence from the written-back `seq.sequence` field,
/// which differs from the phyloxml posterior argmax at gap and unknown positions, so the augur reads
/// are gathered separately from `AncestralOutputMaps`.
#[derive(Debug)]
pub struct AugurOutputMaps {
  /// Reconstructed nucleotide root sequence used as the augur `reference.nuc`.
  pub root_sequence: Seq,
  /// Reconstructed nucleotide sequence per node.
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  /// Parent-edge substitutions per edge (unfiltered; the writer applies the mask filter and pos sort).
  pub edge_subs: BTreeMap<GraphEdgeKey, Vec<Sub>>,
  /// Alignment length (number of sites).
  pub sequence_length: usize,
  /// Ambiguous (unknown) alphabet character used to fill masked positions.
  pub ambiguous_char: AsciiChar,
}

/// Per-node ancestral output as a value: the name and input branch support the output writers read.
#[derive(Debug, Clone, Serialize)]
pub struct AncestralNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

/// Per-edge ancestral output as a value: the branch length the output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}

/// Ancestral reconstruction result as a value.
///
/// The durable per-node and per-edge outputs are reachable directly off the result: `node_sequences`
/// holds the reconstructed sequences captured from the serial reconstruction walk, and the model
/// metadata sits alongside. `graph` carries the tree the output writers still read from.
#[derive(Serialize)]
pub struct AncestralResult {
  #[serde(skip)]
  pub graph: Graph<AncestralGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, AncestralNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
  #[serde(skip)]
  pub aa_node_data: Option<AaNodeData>,
}
