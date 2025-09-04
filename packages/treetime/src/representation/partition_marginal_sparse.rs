use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::GraphEdgeKey;
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::PartitionWithGtrInference;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::graph_sparse::{SparseEdgePartition, SparseNodePartition};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_compressed::PartitionCompressed;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionCompressed for PartitionMarginalSparse {
  fn index(&self) -> usize {
    self.index
  }

  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn length(&self) -> usize {
    self.length
  }

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition> {
    &self.nodes
  }

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &self.edges
  }

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition> {
    &mut self.nodes
  }

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &mut self.edges
  }
}

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionWithGtrInference for PartitionMarginalSparse {
  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn get_seq_composition(&self, node_key: GraphNodeKey) -> &Composition {
    &self.nodes[&node_key].seq.composition
  }

  fn get_edge_substitutions(&self, edge_key: GraphEdgeKey, _graph: &GraphAncestral) -> Vec<Sub> {
    self.edges[&edge_key].subs.clone()
  }
}
