use itertools::Itertools;
use std::collections::{HashSet, VecDeque};

/// Public representation of a node during traversal.
pub trait GraphNode<N> {
  fn index(&self) -> usize;

  fn payload(&self) -> &N;

  fn payload_mut(&mut self) -> &mut N;
}

/// Public representation of a node during forward traversal (from roots to leaves).
///
/// Here we only have access to parent data, because child data is not yet computed.
pub trait GraphNodeForward<'i, N, GN, NI>: GraphNode<N>
where
  N: 'i,
  GN: 'i + GraphNode<N>,
  NI: Iterator<Item = &'i dyn GraphNode<N>>,
{
  fn parents(&self) -> NI;

  fn is_root(&self) -> bool;
}

/// Public representation of a node during reverse traversal (from leaves to root)
///
/// Here we only have access to child data, because parent data is not yet computed.
pub trait GraphNodeReverse<'i, N, GN, NI>: GraphNode<N>
where
  N: 'i,
  GN: 'i + GraphNode<N>,
  NI: Iterator<Item = &'i dyn GraphNode<N>>,
{
  fn children(&self) -> NI;

  fn is_leaf(&self) -> bool;
}

/// Internal representation of a node
#[derive(Clone)]
struct GraphNodeImpl<'g, N> {
  idx: usize,
  parents: Vec<usize>,
  children: Vec<usize>,
  visited: bool,
  payload: N,
  nodes: &'g [GraphNodeImpl<'g, N>],
}

impl<'g, N> GraphNode<N> for GraphNodeImpl<'g, N> {
  fn index(&self) -> usize {
    self.idx
  }

  fn payload(&self) -> &N {
    &self.payload
  }

  fn payload_mut(&mut self) -> &mut N {
    &mut self.payload
  }
}

impl<'g, 'i, N, GN, NI> GraphNodeForward<'i, N, GN, NI> for GraphNodeImpl<'g, N>
where
  N: 'i,
  GN: 'i + GraphNode<N>,
  NI: Iterator<Item = &'i dyn GraphNode<N>>,
{
  fn parents(&self) -> NI {
    self
      .children
      .iter()
      .map(|&idx| &self.nodes[idx] as &'i dyn GraphNode<N>)
  }

  fn is_root(&self) -> bool {
    self.parents.is_empty()
  }
}

impl<'g, 'i, N, GN, NI> GraphNodeReverse<'i, N, GN, NI> for GraphNodeImpl<'g, N>
where
  N: 'i,
  GN: 'i + GraphNode<N>,
  NI: Iterator<Item = &'i dyn GraphNode<N>>,
{
  fn children(&self) -> NI {
    self
      .children
      .iter()
      .map(|&idx| &self.nodes[idx] as &'i dyn GraphNode<N>)
      .into_iter()
  }

  fn is_leaf(&self) -> bool {
    self.children.is_empty()
  }
}

impl<'g, N> GraphNodeImpl<'g, N> {
  pub const fn new(idx: usize, payload: N, nodes: &'g [GraphNodeImpl<'g, N>]) -> Self {
    Self {
      idx,
      parents: vec![],
      children: vec![],
      payload,
      visited: false,
      nodes,
    }
  }

  pub fn is_root(&self) -> bool {
    self.parents.is_empty()
  }

  pub fn is_leaf(&self) -> bool {
    self.children.is_empty()
  }
}

/// Representation of a graph
#[derive(Clone)]
pub struct Graph<'g, N> {
  nodes: Vec<GraphNodeImpl<'g, N>>,
  roots: Vec<usize>,
  leaves: Vec<usize>,
}

impl<'i, 'g, N> Graph<'g, N> {
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
  pub fn insert(&'g mut self, node_payload: N) -> usize {
    let node = GraphNodeImpl::<'g, N>::new(self.nodes.len(), node_payload, self.nodes.as_slice());
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
  pub fn iter_breadth_first_forward<GN, NI, F>(&mut self, visitor_func: F)
  where
    N: 'i,
    GN: 'i + GraphNode<N>,
    NI: Iterator<Item = &'i dyn GraphNode<N>>,
    F: FnMut(&mut dyn GraphNodeForward<'i, N, GN, NI>),
  {
    let mut queue: VecDeque<usize> = self.roots.iter().copied().collect();
    let mut visited = HashSet::<usize>::new();
    while let Some(idx) = queue.pop_front() {
      if !visited.contains(&idx) {
        let node = &mut self.nodes[idx];
        queue.extend(node.children.iter().copied());

        visitor_func(node as &mut dyn GraphNodeForward<'i, N, GN, NI>);

        // visitor_func(&GraphNodeForward {
        //   idx,
        //   parents: node.parents.iter().map(|&idx| &self.nodes[idx].payload).collect_vec(),
        //   payload: &mut node.payload,
        // });
        visited.insert(idx);
      }
    }
  }

  // /// Iterates graph from leaves to roots in breadth-first order
  // ///
  // /// Guarantees that for each visited node, all of it children (recursively) are visited before
  // /// the node itself is visited.
  // ///
  // /// The order of traversal does not take into account the order of insertion of nodes or edges.
  // pub fn iter_breadth_first_reverse(&mut self, visitor_func: impl Fn(&mut dyn GraphNodeReverse<N>)) {
  //   let mut queue: VecDeque<usize> = self.leaves.iter().copied().collect();
  //   let mut visited = HashSet::<usize>::new();
  //   while let Some(idx) = queue.pop_front() {
  //     if !visited.contains(&idx) {
  //       let node = &mut self.nodes[idx];
  //       queue.extend(node.parents.iter().copied());
  //       visitor_func(&mut node);
  //       visited.insert(idx);
  //     }
  //   }
  // }
}
