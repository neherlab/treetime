use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use serde::Serialize;
use treetime_graph::graph::Graph;

pub type GraphTimetree<D = ()> = Graph<NodeTimetree, EdgeTimetree, D>;

#[derive(Debug, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PartitionTimetree {
  Dense(PartitionMarginalDense),
  Sparse(PartitionMarginalSparse),
}
