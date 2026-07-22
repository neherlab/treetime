use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;

#[derive(Serialize)]
pub struct OptimizeGraphData {
  pub gtr: GTR,
  pub model_name: GtrModelName,
  pub sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  pub dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>,
}

impl OptimizeGraphData {
  pub fn new(
    gtr: GTR,
    model_name: GtrModelName,
    sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
    dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>,
  ) -> Self {
    Self {
      gtr,
      model_name,
      sparse_partitions,
      dense_partitions,
    }
  }
}

#[derive(Serialize)]
pub struct OptimizeResult {
  #[serde(skip)]
  pub graph: GraphAncestral<OptimizeGraphData>,
}
