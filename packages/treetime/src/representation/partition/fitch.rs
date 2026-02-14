use crate::alphabet::alphabet::Alphabet;
use crate::representation::partition::traits::PartitionCompressed;
use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Debug)]
pub struct PartitionFitch {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionCompressed for PartitionFitch {
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
