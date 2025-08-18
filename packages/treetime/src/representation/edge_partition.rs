use crate::representation::graph_dense::DenseEdgePartition;
use crate::representation::graph_meta::MetaEdgePartition;
use crate::representation::graph_sparse::SparseEdgePartition;
use enum_extract_macro::EnumExtract;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize, EnumExtract)]
pub enum EdgePartition {
  Dense(DenseEdgePartition),
  Sparse(SparseEdgePartition),
  Meta(MetaEdgePartition),
}

impl EdgePartition {
  pub fn dense() -> Result<Self, Report> {
    Ok(EdgePartition::Dense(DenseEdgePartition::default()))
  }

  pub fn sparse() -> Result<Self, Report> {
    Ok(EdgePartition::Sparse(SparseEdgePartition::default()))
  }

  pub fn meta() -> Result<Self, Report> {
    Ok(EdgePartition::Meta(MetaEdgePartition::new()?))
  }
}
