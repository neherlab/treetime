use std::collections::BTreeMap;

use super::graph_dense::{DenseSeqEdge, DenseSeqNode};
use crate::graph::edge::GraphEdgeKey;
use crate::gtr::gtr::GTR;
use crate::{alphabet::alphabet::Alphabet, graph::node::GraphNodeKey};

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, DenseSeqNode>,
  pub edges: BTreeMap<GraphEdgeKey, DenseSeqEdge>,
}
