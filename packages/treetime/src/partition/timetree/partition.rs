use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;
use treetime_graph::graph::Graph;

pub type GraphTimetree<D = ()> = Graph<NodeTimetree, EdgeTimetree, D>;
pub type PartitionTimetreeRef = Arc<RwLock<PartitionTimetree>>;
pub type PartitionTimetreeAllVec = Vec<PartitionTimetreeRef>;

#[derive(Debug, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PartitionTimetree {
  Dense(PartitionMarginalDense),
  Sparse(PartitionMarginalSparse),
}
