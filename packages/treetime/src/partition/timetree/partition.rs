use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use serde::Serialize;

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PartitionTimetree {
  Dense(DenseReconstruction),
  Sparse(SparseReconstruction),
}
