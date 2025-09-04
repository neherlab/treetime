use super::graph_dense::{DenseEdgePartition, DenseNodePartition};
use crate::graph::edge::GraphEdgeKey;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::PartitionWithGtrInference;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::log_lh::HasLogLh;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
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

impl PartitionWithGtrInference for PartitionMarginalDense {
  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn get_seq_composition(&self, node_key: GraphNodeKey) -> &Composition {
    unimplemented!("Dense node composition is not implemented yet")
  }

  fn get_edge_substitutions(&self, edge_key: GraphEdgeKey, graph: &GraphAncestral) -> Vec<Sub> {
    unimplemented!("Dense node substitutions lookup is not implemented yet")
  }
}
