use crate::graph::core::{CLOSED, OPEN};
use crate::graph::edge::Edge;
use parking_lot::{Mutex, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Weak};

type Outbound<N, E> = RwLock<Vec<Arc<Edge<N, E>>>>;
type Inbound<N, E> = RwLock<Vec<Weak<Edge<N, E>>>>;

/// Represents a node in the graph. Data can be stored in and loaded from the
/// node in a thread safe manner.
#[derive(Debug)]
pub struct Node<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  key: usize,
  pub data: Mutex<N>,
  pub outbound: Outbound<N, E>,
  pub inbound: Inbound<N, E>,
  lock: AtomicBool,
}

impl<N, E> Node<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Create a new node.
  #[inline]
  pub const fn new(key: usize, data: N) -> Node<N, E> {
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
  pub const fn key(&self) -> usize {
    self.key
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

  // /// Find an outbound node and return the corresponding edge if found.
  // #[inline]
  // pub fn find_outbound(&self, target: &Arc<Mutex<Node<N, E>>>) -> Option<Arc<Edge<N, E>>> {
  //   for edge in self.outbound().iter() {
  //     if edge.target().lock().key() == target.key() {
  //       return Some(Arc::clone(edge));
  //     }
  //   }
  //   None
  // }
  //
  // /// Find an inbound node and return the corresponding edge if found.
  // #[inline]
  // pub fn find_inbound(&self, source: &Arc<Mutex<Node<N, E>>>) -> Option<Weak<Edge<N, E>>> {
  //   for edge in self.inbound().iter() {
  //     if edge.upgrade().unwrap().source().lock().key == source.key {
  //       return Some(Weak::clone(edge));
  //     }
  //   }
  //   None
  // }

  /// Get read access to outbound edges of the node.
  #[inline]
  pub fn outbound(&self) -> RwLockReadGuard<Vec<Arc<Edge<N, E>>>> {
    self.outbound.read()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn outbound_mut(&self) -> RwLockWriteGuard<Vec<Arc<Edge<N, E>>>> {
    self.outbound.write()
  }

  /// Get read access to inbound edges of the node.
  #[inline]
  pub fn inbound(&self) -> RwLockReadGuard<Vec<Weak<Edge<N, E>>>> {
    self.inbound.read()
  }

  /// Get read and write access to the outbound edges of the node. Will block other threads.
  #[inline]
  pub fn inbound_mut(&self) -> RwLockWriteGuard<Vec<Weak<Edge<N, E>>>> {
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
  pub fn map_adjacent_dir<F>(&self, user_closure: &F) -> Vec<Weak<Edge<N, E>>>
  where
    N: Clone + Debug + Display + Sync + Send,
    E: Clone + Debug + Display + Sync + Send,
    F: Fn(&Arc<Edge<N, E>>),
  {
    let mut segment: Vec<Weak<Edge<N, E>>> = Vec::new();
    for edge in self.outbound().iter() {
      let target = edge.target();
      if target.read().try_lock() == OPEN {
        target.read().close();
        user_closure(edge);
        segment.push(Arc::downgrade(edge));
      }
    }
    segment
  }
}

impl<N, E> Display for Node<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    let header = format!("{} [label = \"{} : {}\"]", self.key, self.key, self.data.lock());
    write!(fmt, "{}", header)
  }
}
