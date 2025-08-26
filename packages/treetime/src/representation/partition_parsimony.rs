use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::node::GraphNodeKey;
use crate::representation::graph_sparse::{SparseSeqEdge, SparseSeqNode};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct PartitionParsimonyNew {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseSeqNode>,
  pub edges: BTreeMap<GraphEdgeKey, SparseSeqEdge>,
}
