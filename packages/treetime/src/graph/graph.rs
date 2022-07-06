use std::collections::{HashSet, VecDeque};

#[derive(Clone, Debug)]
pub struct GraphNode<N> {
  idx: usize,
  parents: Vec<usize>,
  children: Vec<usize>,
  visited: bool,
  pub payload: N,
}

impl<N> GraphNode<N> {
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

#[derive(Clone, Debug)]
pub struct Graph<N> {
  pub nodes: Vec<GraphNode<N>>,
  pub roots: Vec<usize>,
  pub leaves: Vec<usize>,
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
  pub fn insert(&mut self, node_payload: N) -> usize {
    let node = GraphNode::new(self.nodes.len(), node_payload);
    self.nodes.push(node);
    self.nodes.len() - 1
  }

  /// Adds an edge between two nodes given their indices
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

  /// Iterates graph from roots to leaves in breadth-first order
  pub fn iter_breadth_first(&self, f: impl Fn(&GraphNode<N>)) {
    let mut queue: VecDeque<usize> = self.roots.iter().copied().collect();
    let mut visited = HashSet::<usize>::new();
    while let Some(idx) = queue.pop_front() {
      if !visited.contains(&idx) {
        let node = &self.nodes[idx];
        queue.extend(node.children.iter().copied());
        f(node);
        visited.insert(idx);
      }
    }
  }

  /// Iterates graph from leaves to roots in breadth-first order
  pub fn iter_breadth_first_reverse(&self, f: impl Fn(&GraphNode<N>)) {
    let mut queue: VecDeque<usize> = self.leaves.iter().copied().collect();
    let mut visited = HashSet::<usize>::new();
    while let Some(idx) = queue.pop_front() {
      if !visited.contains(&idx) {
        let node = &self.nodes[idx];
        queue.extend(node.parents.iter().copied());
        f(node);
        visited.insert(idx);
      }
    }
  }
}
