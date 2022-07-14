use crate::graph::node::Node;
use parking_lot::RwLock;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::sync::{Arc, Weak};

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source node's are
/// weak references and can't outlive the nodes they represent.
#[derive(Debug)]
pub struct Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  source: Weak<RwLock<Node<N, E>>>,
  pub target: Weak<RwLock<Node<N, E>>>,
  data: RwLock<E>,
}

impl<N, E> Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Creates a new edge.
  pub fn new(source: Weak<RwLock<Node<N, E>>>, target: Weak<RwLock<Node<N, E>>>, data: E) -> Edge<N, E> {
    Edge {
      source,
      target,
      data: RwLock::new(data),
    }
  }

  /// Edge's source node.
  #[inline]
  pub fn source(&self) -> Arc<RwLock<Node<N, E>>> {
    self.source.upgrade().unwrap()
  }

  /// Edge's target node.
  #[inline]
  pub fn target(&self) -> Arc<RwLock<Node<N, E>>> {
    self.target.upgrade().unwrap()
  }

  /// Load data from the edge.
  #[inline]
  pub fn load(&self) -> E {
    self.data.read().clone()
  }

  /// Store data into the edge.
  #[inline]
  pub fn store(&self, data: E) {
    let mut x = self.data.write();
    *x = data;
  }
}

impl<N, E> Display for Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source().read().key(), self.target().read().key())
  }
}
