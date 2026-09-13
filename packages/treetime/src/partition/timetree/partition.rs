use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use serde::Serialize;

/// A timetree partition as a completed reconstruction bundle: the durable partition inputs together
/// with the node states and per-edge messages/estimates its marginal passes returned. The timetree
/// pipeline threads this and reads branch/optimize data through a transient read view.
#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PartitionTimetree {
  Dense(DenseReconstruction),
  Sparse(SparseReconstruction),
}
