use parking_lot::{Mutex, RwLock, RwLockReadGuard, RwLockWriteGuard};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
  fmt::{Debug, Display, Formatter},
  hash::Hash,
  sync::{
    atomic::{AtomicBool, Ordering},
    Arc, Weak,
  },
};

/// The Traverse enum is used when exploring edges in a graph. It's
/// states signify if an edge should be included in the search, skipped
/// or if the search should stop because a terminating condition has
/// been met such as finding a sink node.
pub enum Traverse {
  Include,
  Skip,
  Finish,
}

/// Represents an empty parameter for either a node or an edge.
#[derive(Clone, Debug)]
pub struct Empty;

impl Display for Empty {
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "_")
  }
}

pub type Frontier<K, N, E> = Vec<Weak<Edge<K, N, E>>>;

pub trait Explorer<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn next_frontier(&self) -> Option<Frontier<K, N, E>>;
  fn prev_frontier(&self) -> Option<Frontier<K, N, E>>;
}

const OPEN: bool = false;
const CLOSED: bool = true;

enum Continue<T> {
  Yes(T),
  No(T),
}

/// Edge representing a connection between two nodes. Relevant data can be
/// stored in the edge atomically. Edge's target and source node's are
/// weak references and can't outlive the nodes they represent.
#[derive(Debug)]
pub struct Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  source: Weak<Node<K, N, E>>,
  target: Weak<Node<K, N, E>>,
  data: Mutex<E>,
  lock: AtomicBool,
}

//=============================================================================

impl<K, N, E> Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Creates a new edge.
  pub fn new(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>, data: E) -> Edge<K, N, E> {
    Edge {
      source: Arc::downgrade(source),
      target: Arc::downgrade(target),
      data: Mutex::new(data),
      lock: AtomicBool::new(OPEN),
    }
  }

  /// Edge's source node.
  #[inline]
  pub fn source(&self) -> Arc<Node<K, N, E>> {
    self.source.upgrade().unwrap()
  }

  /// Edge's target node.
  #[inline]
  pub fn target(&self) -> Arc<Node<K, N, E>> {
    self.target.upgrade().unwrap()
  }

  /// Load data from the edge.
  #[inline]
  pub fn load(&self) -> E {
    self.data.lock().clone()
  }

  /// Store data into the edge.
  #[inline]
  pub fn store(&self, data: E) {
    let mut x = self.data.lock();
    *x = data;
  }

  #[inline]
  fn try_lock(&self) -> bool {
    self.lock.load(Ordering::Relaxed)
  }

  #[inline]
  fn close(&self) {
    self.lock.store(CLOSED, Ordering::Relaxed);
  }

  #[inline]
  fn open(&self) {
    self.lock.store(OPEN, Ordering::Relaxed);
  }
}

//=============================================================================

unsafe impl<K, N, E> Sync for Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
}

impl<K, N, E> Clone for Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn clone(&self) -> Self {
    Edge {
      source: Weak::clone(&self.source),
      target: Weak::clone(&self.target),
      data: Mutex::new(self.data.lock().clone()),
      lock: AtomicBool::new(OPEN),
    }
  }
}

impl<K, N, E> Display for Edge<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn fmt(&self, fmt: &mut Formatter<'_>) -> std::fmt::Result {
    write!(fmt, "{} -> {}", self.source().key(), self.target().key())
  }
}

/// Backtrack edges in an edge list starting from the last edge in the list.
/// Used for example to find the shortest path from the results of a breadth
/// first traversal.
pub fn backtrack_edges<K, N, E>(edges: &Vec<Weak<Edge<K, N, E>>>) -> Vec<Weak<Edge<K, N, E>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  let mut res = Vec::new();
  if edges.is_empty() {
    return res;
  }
  let w = &edges[edges.len() - 1];
  res.push(Weak::clone(w));
  let mut i = 0;
  for edge in edges.iter().rev() {
    let source = res[i].upgrade().unwrap().source();
    if edge.upgrade().unwrap().target() == source {
      res.push(Weak::clone(edge));
      i += 1;
    }
  }
  res.reverse();
  res
}

// Opens all locks in
fn open_locks<K, N, E>(edges: &Vec<Weak<Edge<K, N, E>>>)
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  if edges.is_empty() {
    return;
  }
  edges[0].upgrade().unwrap().source().open();
  for weak in edges.iter() {
    let edge = weak.upgrade().unwrap();
    edge.open();
    edge.target().open();
  }
}

//=============================================================================
// NODE IMPLEMENTATION
//=============================================================================

//=============================================================================
// TYPES

type Outbound<K, N, E> = RwLock<Vec<Arc<Edge<K, N, E>>>>;
type Inbound<K, N, E> = RwLock<Vec<Weak<Edge<K, N, E>>>>;

//=============================================================================
// STRUCT

/// Represents a node in the graph. Data can be stored in and loaded from the
/// node in a thread safe manner.
///
#[derive(Debug)]
pub struct Node<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  key: K,
  data: Mutex<N>,
  outbound: Outbound<K, N, E>,
  inbound: Inbound<K, N, E>,
  lock: AtomicBool,
}

//=============================================================================
// STRUCT IMPLEMENTATIONS

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
      if edge.target() == *target {
        return Some(Arc::clone(edge));
      }
    }
    None
  }

  /// Find an inbound node and return the corresponding edge if found.
  #[inline]
  pub fn find_inbound(&self, source: &Arc<Node<K, N, E>>) -> Option<Weak<Edge<K, N, E>>> {
    for edge in self.inbound().iter() {
      if edge.upgrade().unwrap().source() == *source {
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
  fn close(&self) {
    self.lock.store(CLOSED, Ordering::Relaxed);
  }

  #[inline]
  fn open(&self) {
    self.lock.store(OPEN, Ordering::Relaxed);
  }

  #[inline]
  fn map_adjacent_dir<F>(&self, user_closure: &F) -> Continue<Vec<Weak<Edge<K, N, E>>>>
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

  #[inline]
  fn map_adjacent_undir<F>(&self, user_closure: &F) -> Continue<Vec<Weak<Edge<K, N, E>>>>
  where
    K: Hash + Eq + Clone + Debug + Display + Sync + Send,
    N: Clone + Debug + Display + Sync + Send,
    E: Clone + Debug + Display + Sync + Send,
    F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    let mut segment: Vec<Weak<Edge<K, N, E>>> = Vec::new();
    for edge in self.outbound().iter() {
      if edge.try_lock() == OPEN && edge.target().try_lock() == OPEN {
        edge.close();
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
    for edge in self.inbound().iter() {
      let upgrade = edge.upgrade().unwrap();
      if upgrade.try_lock() == OPEN && upgrade.target().try_lock() == OPEN {
        upgrade.close();
        upgrade.target().close();
        let traversal_state = user_closure(&upgrade);
        match traversal_state {
          Traverse::Include => {
            segment.push(Weak::clone(edge));
          }
          Traverse::Finish => {
            segment.push(Weak::clone(edge));
            return Continue::No(segment);
          }
          Traverse::Skip => {
            upgrade.target().open();
            upgrade.open();
          }
        }
      }
    }
    Continue::Yes(segment)
  }
}

// unsafe impl<K, N, E> Sync for Node<K, N, E>
// where
//   K: Hash + Eq + Clone + Debug + Display + Sync + Send,
//   N: Clone + Debug + Display + Sync + Send,
//   E: Clone + Debug + Display + Sync + Send,
// {
// }

impl<K, N, E> Clone for Node<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn clone(&self) -> Self {
    Node {
      key: self.key.clone(),
      data: Mutex::new(self.data.lock().clone()),
      outbound: Outbound::new(Vec::new()),
      inbound: Inbound::new(Vec::new()),
      lock: AtomicBool::new(OPEN),
    }
  }
}

impl<K, N, E> PartialEq for Node<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn eq(&self, other: &Self) -> bool {
    if self.key == other.key {
      return true;
    }
    false
  }
}

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

#[inline]
fn overlaps<K, N, E>(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  for edge in source.outbound().iter() {
    if edge.target() == *target {
      return true;
    }
  }
  false
}

/// Connect two nodes if no previous connection exists.
pub fn connect<K, N, E>(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>, data: E) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  if !overlaps(source, target) {
    let new_edge = Arc::new(Edge::new(source, target, data));
    source.outbound_mut().push(Arc::clone(&new_edge));
    target.inbound_mut().push(Arc::downgrade(&new_edge));
    return true;
  }
  false
}

/// Disconnect two nodes from each other if they share an edge.
pub fn disconnect<K, N, E>(source: &Arc<Node<K, N, E>>, target: &Arc<Node<K, N, E>>) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  let mut idx: (usize, usize) = (0, 0);
  let mut flag = false;
  for (i, edge) in source.outbound().iter().enumerate() {
    if edge.target() == *target {
      idx.0 = i;
      flag = true;
    }
  }
  for (i, edge) in target.inbound().iter().enumerate() {
    if edge.upgrade().unwrap().source() == *source {
      idx.1 = i;
    }
  }
  if flag {
    source.outbound_mut().remove(idx.0);
    source.inbound_mut().remove(idx.1);
  }
  flag
}

/// # Breadth Traversal
///
/// Conduct a breadth first traversal starting from the source node.
/// User provides an `explorer` closure which determines how nodes and edges
/// are to be interpreted. The closure will return a Traverse enum which has
/// 3 states Include, Skip and Finish. These states determine if we are to
/// "go through" the edge and thus include it in our search. Include will
/// include the edge and continue the search, Skip will indicate that the edge
/// is not to be traversed and Finish will include the edge and finish the
/// algorithm.
///
/// Function will return an `Option<Vec<Weak<Edge<K, N, E>>>>` where a Some value
/// indicates that the traversal was successful ie. a Finish condition was
/// reached. And WeakEdges is a collection of all the traversed edges.
/// The last edge will contain the result that triggered the Finish condition.
/// To get the shortest path for example, we'd backtrack the WeakEdges starting
/// from the last edge which would contain our sink node.
pub fn directed_breadth_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
{
  let mut frontiers: Vec<Weak<Edge<K, N, E>>>;
  let mut bounds: (usize, usize) = (0, 0);
  source.close();
  let initial = source.map_adjacent_dir(&explorer);
  match initial {
    Continue::No(segment) => {
      open_locks(&segment);
      return Some(segment);
    }
    Continue::Yes(segment) => {
      frontiers = segment;
    }
  }
  loop {
    bounds.1 = frontiers.len();
    if bounds.0 == bounds.1 {
      break;
    }
    let current_frontier = &frontiers[bounds.0..bounds.1];
    bounds.0 = bounds.1;
    let mut new_segments = Vec::new();
    for edge in current_frontier.iter() {
      let node = edge.upgrade().unwrap().target();
      let haystack = node.map_adjacent_dir(&explorer);
      match haystack {
        Continue::No(mut segment) => {
          new_segments.append(&mut segment);
          frontiers.append(&mut new_segments);
          open_locks(&frontiers);
          return Some(frontiers);
        }
        Continue::Yes(mut segment) => {
          new_segments.append(&mut segment);
        }
      }
    }
    frontiers.append(&mut new_segments);
  }
  open_locks(&frontiers);
  None
}

pub fn undirected_breadth_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
{
  let mut frontiers: Vec<Weak<Edge<K, N, E>>>;
  let mut bounds: (usize, usize) = (0, 0);
  source.close();
  let initial = source.map_adjacent_undir(&explorer);
  match initial {
    Continue::No(segment) => {
      open_locks(&segment);
      return Some(segment);
    }
    Continue::Yes(segment) => {
      frontiers = segment;
    }
  }
  loop {
    bounds.1 = frontiers.len();
    if bounds.0 == bounds.1 {
      break;
    }
    let current_frontier = &frontiers[bounds.0..bounds.1];
    bounds.0 = bounds.1;
    let mut new_segments = Vec::new();
    for edge in current_frontier.iter() {
      let node = edge.upgrade().unwrap().target();
      let haystack = node.map_adjacent_undir(&explorer);
      match haystack {
        Continue::No(mut segment) => {
          new_segments.append(&mut segment);
          frontiers.append(&mut new_segments);
          open_locks(&frontiers);
          return Some(frontiers);
        }
        Continue::Yes(mut segment) => {
          new_segments.append(&mut segment);
        }
      }
    }
    frontiers.append(&mut new_segments);
  }
  open_locks(&frontiers);
  None
}

/// # Parallel Breadth Traversal
///
/// Conduct a parallel breadth first traversal starting from the source node.
/// User provides an `explorer` closure which determines how nodes and edges
/// are to be interpreted. The closure will return a Traverse enum which has
/// 3 states Include, Skip and Finish. These states determine if we are to
/// "go through" the edge and thus include it in our search. Include will
/// include the edge and continue the search, Skip will indicate that the edge
/// is not to be traversed and Finish will include the edge and finish the
/// algorithm.
///
/// Function will return an `Option<Vec<Weak<Edge<K, N, E>>>>` where a Some value
/// indicates that the traversal was successful ie. a Finish condition was
/// reached. And WeakEdges is a collection of all the traversed edges.
/// The last edge will contain the result that triggered the Finish condition.
/// To get the shortest path for example, we'd backtrack the WeakEdges starting
/// from the last edge which would contain our sink node.
pub fn parallel_directed_breadth_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
{
  let mut frontiers: Vec<Weak<Edge<K, N, E>>>;
  let mut bounds: (usize, usize) = (0, 0);
  let terminate: Arc<AtomicBool> = Arc::new(AtomicBool::new(false));
  source.close();
  match source.map_adjacent_dir(&explorer) {
    Continue::No(segment) => {
      open_locks(&segment);
      return Some(segment);
    }
    Continue::Yes(segment) => {
      frontiers = segment;
    }
  }
  loop {
    bounds.1 = frontiers.len();
    if bounds.0 == bounds.1 {
      break;
    }
    let current_frontier = &frontiers[bounds.0..bounds.1];
    bounds.0 = bounds.1;
    let frontier_segments: Vec<_> = current_frontier
      .into_par_iter()
      .map(|edge| {
        if terminate.load(Ordering::Relaxed) {
          None
        } else {
          let node = edge.upgrade().unwrap().target();
          match node.map_adjacent_dir(&explorer) {
            Continue::No(segment) => {
              terminate.store(true, Ordering::Relaxed);
              Some(segment)
            }
            Continue::Yes(segment) => Some(segment),
          }
        }
      })
      .while_some()
      .collect();
    for mut segment in frontier_segments {
      frontiers.append(&mut segment);
    }
    if terminate.load(Ordering::Relaxed) {
      break;
    }
  }
  open_locks(&frontiers);
  terminate.load(Ordering::Relaxed).then(|| frontiers)
}

pub fn parallel_undirected_breadth_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
{
  let mut frontiers: Vec<Weak<Edge<K, N, E>>>;
  let mut bounds: (usize, usize) = (0, 0);
  let terminate: Arc<AtomicBool> = Arc::new(AtomicBool::new(false));
  source.close();
  match source.map_adjacent_undir(&explorer) {
    Continue::No(segment) => {
      open_locks(&segment);
      return Some(segment);
    }
    Continue::Yes(segment) => {
      frontiers = segment;
    }
  }
  loop {
    bounds.1 = frontiers.len();
    if bounds.0 == bounds.1 {
      break;
    }
    let current_frontier = &frontiers[bounds.0..bounds.1];
    bounds.0 = bounds.1;
    let frontier_segments: Vec<_> = current_frontier
      .into_par_iter()
      .map(|edge| {
        if terminate.load(Ordering::Relaxed) {
          None
        } else {
          let node = edge.upgrade().unwrap().target();
          match node.map_adjacent_undir(&explorer) {
            Continue::No(segment) => {
              terminate.store(true, Ordering::Relaxed);
              Some(segment)
            }
            Continue::Yes(segment) => Some(segment),
          }
        }
      })
      .while_some()
      .collect();
    for mut segment in frontier_segments {
      frontiers.append(&mut segment);
    }
    if terminate.load(Ordering::Relaxed) {
      break;
    }
  }
  open_locks(&frontiers);
  terminate.load(Ordering::Relaxed).then(|| frontiers)
}

//=============================================================================

/// # Depth First Traversal
///
/// Conduct a depth first traversal starting from the source node.
/// User provides an `explorer` closure which determines how nodes and edges
/// are to be interpreted. The closure will return a Traverse enum which has
/// 3 states Include, Skip and Finish. These states determine if we are to
/// "go through" the edge and thus include it in our search. Include will
/// include the edge and continue the search, Skip will indicate that the edge
/// is not to be traversed and Finish will include the edge and finish the
/// algorithm.
///
/// Function will return an `Option<Vec<Weak<Edge<K, N, E>>>>` where a Some value
/// indicates that the traversal was successful ie. a Finish condition was
/// reached. And WeakEdges is a collection of all the traversed edges.
/// The last edge will contain the result that triggered the Finish condition.
/// To get the shortest path for example, we'd backtrack the WeakEdges starting
/// from the last edge which would contain our sink node.
fn directed_depth_traversal_recursion<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  results: &mut Vec<Weak<Edge<K, N, E>>>,
  explorer: F,
) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse,
{
  source.close();
  for edge in source.outbound().iter() {
    if edge.target().try_lock() == OPEN {
      edge.target().close();
      let traverse = explorer(edge);
      match traverse {
        Traverse::Include => {
          results.push(Arc::downgrade(edge));
        }
        Traverse::Finish => {
          results.push(Arc::downgrade(edge));
          return true;
        }
        Traverse::Skip => {
          edge.target().open();
        }
      }
      return directed_depth_traversal_recursion(&edge.target(), results, explorer);
    }
  }
  false
}

pub fn directed_depth_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse,
{
  let mut result = Vec::new();
  let res = directed_depth_traversal_recursion(source, &mut result, explorer);
  open_locks(&result);
  res.then(|| result)
}

fn undirected_depth_traversal_recursion<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  results: &mut Vec<Weak<Edge<K, N, E>>>,
  explorer: F,
) -> bool
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse,
{
  source.close();
  for edge in source.outbound().iter() {
    if edge.target().try_lock() == OPEN {
      edge.target().close();
      let traverse = explorer(edge);
      match traverse {
        Traverse::Include => {
          results.push(Arc::downgrade(edge));
        }
        Traverse::Finish => {
          results.push(Arc::downgrade(edge));
          return true;
        }
        Traverse::Skip => {
          edge.target().open();
        }
      }
      return undirected_depth_traversal_recursion(&edge.target(), results, explorer);
    }
  }
  for edge in source.inbound().iter() {
    let upgrade = edge.upgrade().unwrap();
    if upgrade.target().try_lock() == OPEN {
      upgrade.target().close();
      let traversal_state = explorer(&upgrade);
      match traversal_state {
        Traverse::Include => {
          results.push(Weak::clone(edge));
        }
        Traverse::Finish => {
          results.push(Weak::clone(edge));
          return true;
        }
        Traverse::Skip => {
          upgrade.target().open();
        }
      }
    }
  }
  false
}

pub fn undirected_depth_traversal<K, N, E, F>(
  source: &Arc<Node<K, N, E>>,
  explorer: F,
) -> Option<Vec<Weak<Edge<K, N, E>>>>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
  F: Fn(&Arc<Edge<K, N, E>>) -> Traverse,
{
  let mut result = Vec::new();
  let res = undirected_depth_traversal_recursion(source, &mut result, explorer);
  open_locks(&result);
  res.then(|| result)
}
