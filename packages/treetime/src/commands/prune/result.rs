use crate::gtr::gtr::GTR;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;

#[allow(clippy::manual_non_exhaustive, clippy::partial_pub_fields)] // The private unit field preserves Graph JSON's `data: null` shape.
#[derive(Serialize)]
#[serde(transparent)]
pub struct PruneGraphData {
  marker: (),
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
}

impl PruneGraphData {
  pub fn new(gtr: Option<GTR>, partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>) -> Self {
    Self {
      marker: (),
      gtr,
      partitions,
    }
  }
}

#[derive(Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: GraphAncestral<PruneGraphData>,
}
