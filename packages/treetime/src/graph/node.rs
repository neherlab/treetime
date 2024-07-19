use crate::graph::edge::GraphEdgeKey;
use derive_more::Display;
use eyre::Report;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// Defines how to read and write node name
pub trait Named {
  fn name(&self) -> &str;
  fn set_name(&mut self, name: impl AsRef<str>);
}

/// Defines how to construct node when reading from Newick and Nexus files
pub trait NodeFromNwk: Sized {
  fn from_nwk(name: impl AsRef<str>, comments: &BTreeMap<String, String>) -> Result<Self, Report>;
}

/// Defines how to display node information when writing to Newick and Nexus files
pub trait NodeToNwk {
  fn nwk_name(&self) -> Option<impl AsRef<str>>;

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}

/// Defines how to display node information when writing to GraphViz (.dot) file
pub trait NodeToGraphviz {
  // Defines how to display label (name) of the node in GraphViz (.dot) file
  fn to_graphviz_label(&self) -> String;

  // Defines how to display additional attributes of the node in GraphViz (.dot) file
  fn to_graphviz_attributes(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}

pub trait GraphNode: Clone + Debug + Sync + Send {}

#[derive(Copy, Clone, Debug, Display, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct GraphNodeKey(pub usize);

impl GraphNodeKey {
  #[inline]
  pub const fn as_usize(self) -> usize {
    self.0
  }
}

/// Internal representation of a node in a graph
#[derive(Debug, Serialize, Deserialize)]
pub struct Node<N: GraphNode> {
  key: GraphNodeKey,
  data: Arc<RwLock<N>>,
  outbound_edges: Vec<GraphEdgeKey>,
  inbound_edges: Vec<GraphEdgeKey>,
  is_visited: AtomicBool,
}

impl<N> PartialEq<Self> for Node<N>
where
  N: GraphNode,
{
  fn eq(&self, other: &Self) -> bool {
    self.key == other.key
  }
}

impl<N> Node<N>
where
  N: GraphNode,
{
  /// Create a new node.
  #[inline]
  pub fn new(key: GraphNodeKey, data: N) -> Node<N> {
    Self {
      key,
      data: Arc::new(RwLock::new(data)),
      outbound_edges: Vec::new(),
      inbound_edges: Vec::new(),
      is_visited: AtomicBool::new(false),
    }
  }

  #[inline]
  pub fn payload(&self) -> Arc<RwLock<N>> {
    Arc::clone(&self.data)
  }

  /// Get node key.
  #[inline]
  pub const fn key(&self) -> GraphNodeKey {
    self.key
  }

  /// Get node degree i.e. number of outbound edges.
  #[inline]
  pub fn degree(&self) -> usize {
    self.outbound().len()
  }

  /// Check if node is a leaf node, i.e. has no outbound edges.
  #[inline]
  pub fn is_leaf(&self) -> bool {
    self.outbound().is_empty()
  }

  /// Check if node is a root node, i.e. has no inbound edges.
  #[inline]
  pub fn is_root(&self) -> bool {
    self.inbound().is_empty()
  }

  /// Get read access to outbound edges of the node.
  #[inline]
  pub fn outbound(&self) -> &[GraphEdgeKey] {
    self.outbound_edges.as_slice()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn outbound_mut(&mut self) -> &mut Vec<GraphEdgeKey> {
    &mut self.outbound_edges
  }

  /// Get read access to inbound edges of the node.
  #[inline]
  pub fn inbound(&self) -> &[GraphEdgeKey] {
    self.inbound_edges.as_slice()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn inbound_mut(&mut self) -> &mut Vec<GraphEdgeKey> {
    &mut self.inbound_edges
  }

  #[inline]
  pub fn is_visited(&self) -> bool {
    self.is_visited.load(Ordering::Relaxed)
  }

  #[inline]
  pub fn mark_as_visited(&self) {
    self.is_visited.store(true, Ordering::Relaxed);
  }

  #[inline]
  pub fn mark_as_not_visited(&self) {
    self.is_visited.store(false, Ordering::Relaxed);
  }
}
