use crate::gtr::gtr::GTR;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;

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
  pub seq: Option<Arc<RwLock<PartitionMarginalSparse>>>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
}
