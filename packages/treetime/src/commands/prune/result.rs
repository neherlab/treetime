use crate::gtr::gtr::GTR;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
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

#[derive(Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: GraphAncestral<PruneGraphData>,
}
