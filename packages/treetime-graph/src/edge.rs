use crate::graph::Graph;
use crate::node::GraphNodeKey;
use derive_more::Display;
use getset::{CopyGetters, Getters, MutGetters, Setters};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct GraphEdgeKey(pub usize);

impl GraphEdgeKey {
  pub const fn as_usize(self) -> usize {
    self.0
  }

  pub fn invalid() -> Self {
    Self(usize::MAX)
  }
}

/// Edge representing a connection between two nodes. Edge's target and source
/// nodes are weak references and can't outlive the nodes they represent.
#[derive(Clone, Debug, Serialize, Deserialize, Getters, CopyGetters, MutGetters, Setters)]
pub struct Edge {
  #[getset(get_copy = "pub", get_mut = "pub", set = "pub")]
  key: GraphEdgeKey,

  #[getset(get_copy = "pub", get_mut = "pub", set = "pub")]
  source: GraphNodeKey,

  #[getset(get_copy = "pub", get_mut = "pub", set = "pub")]
  target: GraphNodeKey,
}

impl Edge {
  /// Creates a new edge.
  pub fn new(key: GraphEdgeKey, source: GraphNodeKey, target: GraphNodeKey) -> Edge {
    Edge { key, source, target }
  }
}

/// Invert direction of an edge.
pub fn invert_edge(graph: &mut Graph, edge_key: GraphEdgeKey) {
  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge is not attached to this graph");
    (edge.source(), edge.target())
  };

  // Move this edge from outbound edges to inbound edges of the source node
  {
    let source = graph.get_node_mut(source_key).expect("Edge source node must exist");
    source.outbound_mut().retain(|edge| *edge != edge_key);
    source.inbound_mut().push(edge_key);
  }

  // Move this edge from inbound edges to outbound edges of the target node
  {
    let target = graph.get_node_mut(target_key).expect("Edge target node must exist");
    target.inbound_mut().retain(|edge| *edge != edge_key);
    target.outbound_mut().push(edge_key);
  }

  // Swap source and target nodes inside the edge itself
  {
    let edge = graph.get_edge_mut(edge_key).expect("Edge must exist");
    edge.set_source(target_key);
    edge.set_target(source_key);
  }
}
