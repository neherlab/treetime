use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use serde::Serialize;

#[derive(Debug, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PartitionTimetree {
  Dense(PartitionMarginalDense),
  Sparse(PartitionMarginalSparse),
}
