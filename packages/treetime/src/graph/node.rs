use crate::graph::core::{Continue, Inbound, Outbound, Traverse, CLOSED, OPEN};
use crate::graph::edge::Edge;
use parking_lot::{Mutex, RwLockReadGuard, RwLockWriteGuard};
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Weak};

/// Represents a node in the graph. Data can be stored in and loaded from the
/// node in a thread safe manner.
#[derive(Debug)]
pub struct Node<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub key: K,
  pub data: Mutex<N>,
  pub outbound: Outbound<K, N, E>,
  pub inbound: Inbound<K, N, E>,
  pub lock: AtomicBool,
}

impl<K, N, E> Node<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Create a new node.
  #[inline]
  pub const fn new(key: K, data: N) -> Node<K, N, E> {
    Self {
      key,
      data: Mutex::new(data),
      outbound: Outbound::new(Vec::new()),
      inbound: Inbound::new(Vec::new()),
      lock: AtomicBool::new(OPEN),
    }
  }

  /// Load data from the node.
  #[inline]
  pub fn load(&self) -> N {
    self.data.lock().clone()
  }

  /// Store data to the node.
  #[inline]
  pub fn store(&self, data: N) {
    *self.data.lock() = data;
  }

  /// Get node key.
  #[inline]
  pub const fn key(&self) -> &K {
    &self.key
  }

  /// Get node degree ie. amount of outbound edges.
  #[inline]
  pub fn degree(&self) -> usize {
    self.outbound().len()
  }

  /// Check if node is a leaf node ie. has no outbound edges.
  #[inline]
  pub fn is_leaf(&self) -> bool {
    self.outbound().len() == 0
  }

  /// Find an outbound node and return the corresponding edge if found.
  #[inline]
  pub fn find_outbound(&self, target: &Arc<Node<K, N, E>>) -> Option<Arc<Edge<K, N, E>>> {
    for edge in self.outbound().iter() {
      if edge.target().key() == target.key() {
        return Some(Arc::clone(edge));
      }
    }
    None
  }

  /// Find an inbound node and return the corresponding edge if found.
  #[inline]
  pub fn find_inbound(&self, source: &Arc<Node<K, N, E>>) -> Option<Weak<Edge<K, N, E>>> {
    for edge in self.inbound().iter() {
      if edge.upgrade().unwrap().source().key == source.key {
        return Some(Weak::clone(edge));
      }
    }
    None
  }

  /// Get read access to outbound edges of the node.
  #[inline]
  pub fn outbound(&self) -> RwLockReadGuard<Vec<Arc<Edge<K, N, E>>>> {
    self.outbound.read()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn outbound_mut(&self) -> RwLockWriteGuard<Vec<Arc<Edge<K, N, E>>>> {
    self.outbound.write()
  }

  /// Get read access to inbound edges of the node.
  #[inline]
  pub fn inbound(&self) -> RwLockReadGuard<Vec<Weak<Edge<K, N, E>>>> {
    self.inbound.read()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn inbound_mut(&self) -> RwLockWriteGuard<Vec<Weak<Edge<K, N, E>>>> {
    self.inbound.write()
  }

  #[inline]
  fn try_lock(&self) -> bool {
    self.lock.load(Ordering::Relaxed)
  }

  #[inline]
  pub fn close(&self) {
    self.lock.store(CLOSED, Ordering::Relaxed);
  }

  #[inline]
  pub fn open(&self) {
    self.lock.store(OPEN, Ordering::Relaxed);
  }

  #[inline]
  pub fn map_adjacent_dir<F>(&self, user_closure: &F) -> Continue<Vec<Weak<Edge<K, N, E>>>>
  where
    K: Hash + Eq + Clone + Debug + Display + Sync + Send,
    N: Clone + Debug + Display + Sync + Send,
    E: Clone + Debug + Display + Sync + Send,
    F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    let mut segment: Vec<Weak<Edge<K, N, E>>> = Vec::new();
    for edge in self.outbound().iter() {
      if edge.target().try_lock() == OPEN {
        edge.target().close();
        let traversal_state = user_closure(edge);
        match traversal_state {
          Traverse::Include => {
            segment.push(Arc::downgrade(edge));
          }
          Traverse::Finish => {
            segment.push(Arc::downgrade(edge));
            return Continue::No(segment);
          }
          Traverse::Skip => {
            edge.target().open();
          }
        }
      }
    }
    Continue::Yes(segment)
  }
}

// impl<K, N, E> PartialEq for Node<K, N, E>
// where
//   K: Hash + Eq + Clone + Debug + Display + Sync + Send,
//   N: Clone + Debug + Display + Sync + Send,
//   E: Clone + Debug + Display + Sync + Send,
// {
//   fn eq(&self, other: &Self) -> bool {
//     if self.key == other.key {
//       return true;
//     }
//     false
//   }
// }

impl<K, N, E> Display for Node<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    let header = format!("{} [label = \"{} : {}\"]", self.key, self.key, self.data.lock());
    write!(fmt, "{}", header)
  }
}
