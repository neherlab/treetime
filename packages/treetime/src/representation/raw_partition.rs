use crate::alphabet::alphabet::Alphabet;
use crate::representation::graph_dense::DenseNodePartition;
use crate::representation::graph_meta::MetaNodePartition;
use crate::representation::graph_sparse::SparseNodePartition;
use crate::representation::partitions_likelihood::PartitionLikelihoodWithAln;
use crate::representation::partitions_parsimony::PartitionParsimonyWithAln;
use crate::representation::seq::Seq;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum RawPartition {
  Parsimony(PartitionParsimonyWithAln),
  Likelihood(PartitionLikelihoodWithAln),
}
