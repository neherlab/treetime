use std::collections::VecDeque;

#[derive(Clone, Debug)]
pub struct GraphNode<N> {
  idx: usize,
  parents: Vec<usize>,
  children: Vec<usize>,
  pub payload: N,
}

impl<N> GraphNode<N> {
  pub const fn new(idx: usize, payload: N) -> Self {
    Self {
      idx,
      parents: vec![],
      children: vec![],
      payload,
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
}

impl<N> Graph<N> {
  /// Creates an empty graph
  pub const fn new() -> Self {
    Self { nodes: vec![] }
  }

  /// Returns reference to root node
  pub fn root(&self) -> &GraphNode<N> {
    &self.nodes[0]
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

  /// Iterates graph from root to leaves in breadth-first order
  pub fn iter_breadth_first(&self, f: impl Fn(&GraphNode<N>)) {
    let root = self.root();
    let mut queue = VecDeque::from(vec![root.idx]);

    while let Some(idx) = queue.pop_front() {
      let node = &self.nodes[idx];
      queue.extend(node.children.iter().copied());
      f(node);
    }
  }
}
