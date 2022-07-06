use itertools::Itertools;
use std::collections::{HashSet, VecDeque};

/// Internal representation of a node
#[derive(Clone, Debug)]
struct GraphNodeImpl<N> {
  idx: usize,
  parents: Vec<usize>,
  children: Vec<usize>,
  visited: bool,
  payload: N,
}

impl<N> GraphNodeImpl<N> {
  pub const fn new(idx: usize, payload: N) -> Self {
    Self {
      idx,
      parents: vec![],
      children: vec![],
      payload,
      visited: false,
    }
  }

  pub fn is_root(&self) -> bool {
    self.parents.is_empty()
  }

  pub fn is_leaf(&self) -> bool {
    self.children.is_empty()
  }
}

/// Public representation of a node during forward traversal (from roots to leaves).
///
/// Here we only have access to parent data, because child data is not yet computed.
#[derive(Clone, Debug)]
pub struct GraphNodeForward<'node, N> {
  pub idx: usize,
  pub parents: Vec<&'node N>,
  pub payload: &'node N,
}

impl<'node, N> GraphNodeForward<'node, N> {
  pub fn is_root(&self) -> bool {
    self.parents.is_empty()
  }
}

/// Public representation of a node during reverse traversal
///
/// Here we only have access to child data, because parent data is not yet computed.
#[derive(Clone, Debug)]
pub struct GraphNodeReverse<'node, N> {
  pub idx: usize,
  pub children: Vec<&'node N>,
  pub payload: &'node N,
}

impl<'node, N> GraphNodeReverse<'node, N> {
  pub fn is_leaf(&self) -> bool {
    self.children.is_empty()
  }
}

/// Representation of a graph
#[derive(Clone, Debug)]
pub struct Graph<N> {
  nodes: Vec<GraphNodeImpl<N>>,
  roots: Vec<usize>,
  leaves: Vec<usize>,
}

impl<N> Graph<N> {
  /// Creates an empty graph
  pub const fn new() -> Self {
    Self {
      nodes: vec![],
      roots: vec![],
      leaves: vec![],
    }
  }

  /// Adds a node, given its payload. Returns index of the added node.
  ///
  /// The order of insertion can be arbitrary and is not taken into account in any of the graph
  /// operations.
  pub fn insert(&mut self, node_payload: N) -> usize {
    let node = GraphNodeImpl::new(self.nodes.len(), node_payload);
    self.nodes.push(node);
    self.nodes.len() - 1
  }

  /// Adds an edge between two nodes given their indices
  ///
  /// The edge is directed from the first node to the second node, i.e. the first node is considered
  /// a parent and the second node is considered a child.
  pub fn connect(&mut self, from_idx: usize, to_idx: usize) {
    let from = &mut self.nodes[from_idx];
    from.children.push(to_idx);

    let to = &mut self.nodes[to_idx];
    to.parents.push(from_idx);
  }

  /// Finalizes the graph initialization
  pub fn build(&mut self) -> &mut Self {
    for node in &self.nodes {
      if node.is_leaf() {
        self.leaves.push(node.idx);
      }
      if node.is_root() {
        self.roots.push(node.idx);
      }
    }
    self
  }

  /// Iterates graph from roots to leaves in breadth-first order.
  ///
  /// Guarantees that for each visited node, all of it parents (recursively) are visited before
  /// the node itself is visited.
  ///
  /// The order of traversal does not take into account the order of insertion of nodes or edges.
  pub fn iter_breadth_first_forward(&self, visitor_func: impl Fn(&GraphNodeForward<N>)) {
    let mut queue: VecDeque<usize> = self.roots.iter().copied().collect();
    let mut visited = HashSet::<usize>::new();
    while let Some(idx) = queue.pop_front() {
      if !visited.contains(&idx) {
        let node = &self.nodes[idx];
        queue.extend(node.children.iter().copied());
        visitor_func(&GraphNodeForward {
          idx,
          parents: node.parents.iter().map(|&idx| &self.nodes[idx].payload).collect_vec(),
          payload: &node.payload,
        });
        visited.insert(idx);
      }
    }
  }

  /// Iterates graph from leaves to roots in breadth-first order
  ///
  /// Guarantees that for each visited node, all of it children (recursively) are visited before
  /// the node itself is visited.
  ///
  /// The order of traversal does not take into account the order of insertion of nodes or edges.
  pub fn iter_breadth_first_reverse(&self, visitor_func: impl Fn(&GraphNodeReverse<N>)) {
    let mut queue: VecDeque<usize> = self.leaves.iter().copied().collect();
    let mut visited = HashSet::<usize>::new();
    while let Some(idx) = queue.pop_front() {
      if !visited.contains(&idx) {
        let node = &self.nodes[idx];
        queue.extend(node.parents.iter().copied());
        visitor_func(&GraphNodeReverse {
          idx,
          children: node.children.iter().map(|&idx| &self.nodes[idx].payload).collect_vec(),
          payload: &node.payload,
        });
        visited.insert(idx);
      }
    }
  }
}
