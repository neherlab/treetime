use crate::graph::graph::Weighted;
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
  E: Clone + Debug + Display + Sync + Send + Weighted,
{
  source: Weak<RwLock<Node<N, E>>>,
  target: Weak<RwLock<Node<N, E>>>,
  data: Arc<RwLock<E>>,
}

impl<N, E> Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send + Weighted,
{
  /// Creates a new edge.
  pub fn new(source: Weak<RwLock<Node<N, E>>>, target: Weak<RwLock<Node<N, E>>>, data: E) -> Edge<N, E> {
    Edge {
      source,
      target,
      data: Arc::new(RwLock::new(data)),
    }
  }

  #[inline]
  pub fn source(&self) -> Arc<RwLock<Node<N, E>>> {
    self.source.upgrade().unwrap()
  }

  #[inline]
  pub fn target(&self) -> Arc<RwLock<Node<N, E>>> {
    self.target.upgrade().unwrap()
  }

  #[inline]
  pub fn payload(&self) -> Arc<RwLock<E>> {
    Arc::clone(&self.data)
  }
}

impl<N, E> Display for Edge<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send + Weighted,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source().read().key(), self.target().read().key())
  }
}
