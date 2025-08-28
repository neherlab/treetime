use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::GTR;
use crate::representation::graph_sparse::{SparseSeqEdge, SparseSeqNode};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseSeqNode>,
  pub edges: BTreeMap<GraphEdgeKey, SparseSeqEdge>,
}
