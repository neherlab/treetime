use crate::ancestral::pipeline::AncestralPartition;
use crate::commands::ancestral::aa_node_data::AaNodeData;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::payload::ancestral::GraphAncestral;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

#[derive(Serialize)]
pub struct AncestralGraphData {
  pub partition: Option<AncestralPartition>,
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  pub mask: Vec<bool>,
  pub aa_node_data: Option<AaNodeData>,
}

impl AncestralGraphData {
  pub fn new(
    partition: Option<AncestralPartition>,
    gtr: Option<GTR>,
    model_name: GtrModelName,
    mask: Vec<bool>,
    aa_node_data: Option<AaNodeData>,
  ) -> Self {
    Self {
      partition,
      gtr,
      model_name,
      mask,
      aa_node_data,
    }
  }
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
/// The durable per-node and per-edge outputs are reachable directly off the result: `seq` holds the
/// sequence store (dense, sparse, and parsimony kept as separate representations), `node_sequences`
/// holds the reconstructed sequences captured from the serial reconstruction walk, and the model
/// metadata sits alongside. `graph` carries the tree the output writers still read from; the writers
/// move onto the result value in a later step, after which `graph` and the partition-in-graph go
/// away.
#[derive(Serialize)]
pub struct AncestralResult {
  #[serde(skip)]
  pub graph: GraphAncestral<AncestralGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, AncestralNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub seq: Option<AncestralPartition>,
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
