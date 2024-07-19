use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey};
use derive_more::Display;
use eyre::Report;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::mem::swap;
use std::sync::Arc;

/// Defines how to read and write edge weight
pub trait Weighted {
  fn weight(&self) -> f64;
}

/// Defines how to construct edge when reading from Newick and Nexus files
pub trait EdgeFromNwk: Sized {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report>;
}

/// Defines how to display edge information when writing to Newick and Nexus files
pub trait EdgeToNwk {
  fn nwk_weight(&self) -> Option<f64>;
}

/// Defines how to display edge information when writing to GraphViz (.dot) file
pub trait EdgeToGraphViz {
  // Defines how to display label (name) of the edge in GraphViz (.dot) file
  fn to_graphviz_label(&self) -> String;

  // Defines how to assign weight of the edge in GraphViz (.dot) file
  fn to_graphviz_weight(&self) -> f64;

  // Defines how to display additional attributes of the edge in GraphViz (.dot) file
  fn to_graphviz_attributes(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}

pub trait GraphEdge: Clone + Debug + Sync + Send {}

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct GraphEdgeKey(pub usize);

impl GraphEdgeKey {
  pub const fn as_usize(self) -> usize {
    self.0
  }
}

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source nodes are
/// weak references and can't outlive the nodes they represent.
#[derive(Clone, Debug, Serialize, Deserialize)]
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
  pub const fn key(&self) -> GraphEdgeKey {
    self.key
  }
  pub const fn source(&self) -> GraphNodeKey {
    self.source
  }
  pub const fn target(&self) -> GraphNodeKey {
    self.target
  }
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::examples::tests::get_example_tree;
  use eyre::Report;
  use rstest::rstest;

  #[rstest]
  fn edge_inverts() -> Result<(), Report> {
    let graph = get_example_tree()?;

    let edge = graph.get_edge(GraphEdgeKey(3)).unwrap();
    let input_source = edge.read().source();
    let input_target = edge.read().target();

    invert_edge(&graph, &edge);

    let edge = graph.get_edge(GraphEdgeKey(3)).unwrap();
    let output_source = edge.read().source();
    let output_target = edge.read().target();

    assert_eq!(input_source, output_target);
    assert_eq!(input_target, output_source);

    Ok(())
  }
}
