use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey};
use derive_more::Display;
use parking_lot::RwLock;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::mem::swap;
use std::sync::Arc;

pub trait Weighted {
  fn weight(&self) -> f64;
}

pub trait GraphEdge: Clone + Debug + Display + Sync + Send + Weighted {
  fn new(weight: f64) -> Self;
}

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct GraphEdgeKey(pub usize);

impl GraphEdgeKey {
  #[inline]
  pub const fn as_usize(self) -> usize {
    self.0
  }
}

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source node's are
/// weak references and can't outlive the nodes they represent.
#[derive(Debug)]
pub struct Edge<E: GraphEdge> {
  key: GraphEdgeKey,
  source: GraphNodeKey,
  target: GraphNodeKey,
  data: Arc<RwLock<E>>,
}

impl<E: GraphEdge> Edge<E> {
  /// Creates a new edge.
  pub fn new(key: GraphEdgeKey, source: GraphNodeKey, target: GraphNodeKey, data: E) -> Edge<E> {
    Edge {
      key,
      source,
      target,
      data: Arc::new(RwLock::new(data)),
    }
  }

  #[inline]
  pub const fn key(&self) -> GraphEdgeKey {
    self.key
  }

  #[inline]
  pub const fn source(&self) -> GraphNodeKey {
    self.source
  }

  #[inline]
  pub const fn target(&self) -> GraphNodeKey {
    self.target
  }

  #[inline]
  pub fn payload(&self) -> Arc<RwLock<E>> {
    Arc::clone(&self.data)
  }
}

/// Invert direction of an edge.
pub fn invert_edge<N: GraphNode, E: GraphEdge>(graph: &Graph<N, E>, edge: &Arc<RwLock<Edge<E>>>) {
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
    let edge: &mut Edge<E> = &mut edge.write();
    swap(&mut edge.source, &mut edge.target);
  }
}

impl<E: GraphEdge> Display for Edge<E> {
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source(), self.target())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::examples::get_example_tree;
  use eyre::Report;
  use rstest::rstest;

  #[rstest]
  fn edge_inverts() -> Result<(), Report> {
    let mut graph = get_example_tree()?;

    let edge = graph.get_edge(GraphEdgeKey(3)).unwrap();
    let input_source = edge.read().source();
    let input_target = edge.read().target();

    invert_edge(&mut graph, &edge);

    let edge = graph.get_edge(GraphEdgeKey(3)).unwrap();
    let output_source = edge.read().source();
    let output_target = edge.read().target();

    assert_eq!(input_source, output_target);
    assert_eq!(input_target, output_source);

    Ok(())
  }
}
