use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;

#[allow(clippy::manual_non_exhaustive, clippy::partial_pub_fields)] // The private unit field preserves Graph JSON's `data: null` shape.
#[derive(Serialize)]
#[serde(transparent)]
pub struct OptimizeGraphData {
  marker: (),
  #[serde(skip)]
  pub gtr: GTR,
  #[serde(skip)]
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  #[serde(skip)]
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
      marker: (),
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
