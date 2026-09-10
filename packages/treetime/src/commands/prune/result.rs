use crate::gtr::gtr::GTR;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use serde::Serialize;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

#[derive(Serialize)]
pub struct PruneGraphData {
  pub gtr: Option<GTR>,
  pub partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
}

impl PruneGraphData {
  pub fn new(gtr: Option<GTR>, partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>) -> Self {
    Self { gtr, partitions }
  }
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
/// The durable per-edge outputs are reachable directly off the result: `seq` holds the sparse
/// sequence partition carrying the re-oriented substitutions and indels that the serial topology
/// edits leave on the surviving edges (ids stay stable; a removed node or edge leaves a gap and
/// nothing is renumbered), and `gtr` holds the fitted model. `graph` carries the tree the output
/// writers still read from; the writers move onto the result value in a later step, after which
/// `graph` and the partition-in-graph go away.
#[derive(Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: GraphAncestral<PruneGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, PruneNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub seq: Option<Arc<RwLock<PartitionMarginalSparse>>>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
}
