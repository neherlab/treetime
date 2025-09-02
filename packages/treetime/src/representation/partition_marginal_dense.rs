use super::graph_dense::{DenseEdgePartition, DenseNodePartition};
use crate::graph::edge::GraphEdgeKey;
use crate::gtr::gtr::GTR;
use crate::representation::log_lh::HasLogLh;
use crate::{alphabet::alphabet::Alphabet, graph::node::GraphNodeKey};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, DenseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, DenseEdgePartition>,
}

impl HasLogLh for PartitionMarginalDense {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}
