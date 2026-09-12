use crate::graph::Graph;
use crate::node::GraphNodeKey;
use derive_more::Display;
use getset::{CopyGetters, Getters, MutGetters, Setters};
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;
use std::mem::swap;
use std::sync::Arc;

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
pub fn invert_edge(graph: &mut Graph, edge: &Arc<RwLock<Edge>>) {
  let (this_edge_key, source, target) = {
    let edge = edge.read();

    let this_edge_key = edge.key();

    let source = graph
      .get_node(edge.source())
      .expect("Edge is not attached to this graph");

    let target = graph.get_node(edge.target()).expect("edge must have a target node");

    (this_edge_key, source, target)
  };

  // Move this edge from outbound edges to inbound edges of the source node
  source.write().outbound_mut().retain(|edge| *edge != this_edge_key);
  source.write().inbound_mut().push(this_edge_key);

  // Move this edge from inbound edges to outbound edges of the target node
  target.write().inbound_mut().retain(|edge| *edge != this_edge_key);
  target.write().outbound_mut().push(this_edge_key);

  // Swap source and target nodes inside the edge itself
  {
    let edge: &mut Edge = &mut edge.write();
    swap(&mut edge.source, &mut edge.target);
  }
}
