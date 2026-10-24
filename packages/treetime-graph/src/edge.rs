use crate::graph::Graph;
use crate::node::GraphNodeKey;
use derive_more::Display;
use deser::Serialize;
use getset::{CopyGetters, Getters, MutGetters, Setters};
use std::fmt::Debug;
use std::hash::Hash;

#[derive(Clone, Debug, Getters, CopyGetters, MutGetters, Setters, Serialize)]
pub struct Edge {
  #[getset(get_copy = "pub", get_mut = "pub", set = "pub")]
  key: GraphEdgeKey,

  #[getset(get_copy = "pub", get_mut = "pub", set = "pub")]
  source: GraphNodeKey,

  #[getset(get_copy = "pub", get_mut = "pub", set = "pub")]
  target: GraphNodeKey,
}

impl Edge {
  pub(crate) fn new(key: GraphEdgeKey, source: GraphNodeKey, target: GraphNodeKey) -> Edge {
    Edge { key, source, target }
  }
}

pub fn invert_edge(graph: &mut Graph, edge_key: GraphEdgeKey) {
  let (source_key, target_key) = {
    let edge = &graph[edge_key];
    (edge.source(), edge.target())
  };

  {
    let source = graph.node_mut(source_key);
    source.outbound_mut().retain(|edge| *edge != edge_key);
    source.inbound_mut().push(edge_key);
  }

  {
    let target = graph.node_mut(target_key);
    target.inbound_mut().retain(|edge| *edge != edge_key);
    target.outbound_mut().push(edge_key);
  }

  {
    let edge = graph.edge_mut(edge_key);
    edge.set_source(target_key);
    edge.set_target(source_key);
  }
}

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize)]
pub struct GraphEdgeKey(pub usize);

impl GraphEdgeKey {
  pub const fn as_usize(self) -> usize {
    self.0
  }
}
