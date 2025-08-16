use crate::representation::graph_dense::DenseEdgePartition;
use crate::representation::graph_meta::MetaEdgePartition;
use crate::representation::graph_sparse::SparseEdgePartition;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum EdgePartition {
  Dense(DenseEdgePartition),
  Sparse(SparseEdgePartition),
  Meta(MetaEdgePartition),
}
