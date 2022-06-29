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

  /// Returns iterator to breadth-first traversal
  pub fn iter_breadth_first<'a>(&'a self) -> GraphBreadthFirstIterator<'a, N> {
    GraphBreadthFirstIterator::<'a, N>::new(self)
  }
}

#[derive(Clone, Debug)]
pub struct GraphBreadthFirstIterator<'a, T> {
  graph: &'a Graph<T>,
  queue: VecDeque<usize>,
}

impl<'a, T> GraphBreadthFirstIterator<'a, T> {
  #[inline]
  pub fn new(graph: &'a Graph<T>) -> Self {
    let root = graph.root();
    Self {
      graph,
      queue: VecDeque::from(vec![root.idx]),
    }
  }
}

impl<'a, T> Iterator for GraphBreadthFirstIterator<'a, T> {
  type Item = &'a GraphNode<T>;

  #[inline]
  fn next(&mut self) -> Option<Self::Item> {
    self.queue.pop_front().map(|idx| {
      let node = &self.graph.nodes[idx];
      self.queue.extend(node.children.iter().copied());
      node
    })
  }
}
