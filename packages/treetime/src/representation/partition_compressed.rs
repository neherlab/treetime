use crate::alphabet::alphabet::Alphabet;
use crate::representation::graph_sparse::{SparseEdgePartition, SparseNodePartition};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

pub trait PartitionCompressed: Sync + Send {
  fn index(&self) -> usize;

  fn alphabet(&self) -> &Alphabet;

  fn length(&self) -> usize;

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition>;

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition>;

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition>;

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition>;

  fn node(&self, key: &GraphNodeKey) -> &SparseNodePartition {
    self.nodes().get(key).expect("Node not found")
  }

  fn edge(&self, key: &GraphEdgeKey) -> &SparseEdgePartition {
    self.edges().get(key).expect("Edge not found")
  }

  fn node_mut(&mut self, key: &GraphNodeKey) -> &mut SparseNodePartition {
    self.nodes_mut().get_mut(key).expect("Node not found")
  }

  fn edge_mut(&mut self, key: &GraphEdgeKey) -> &mut SparseEdgePartition {
    self.edges_mut().get_mut(key).expect("Edge not found")
  }
}
