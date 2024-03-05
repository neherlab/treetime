use crate::graph::breadth_first::{
  directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward, GraphTraversalContinuation,
};
use crate::graph::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::node::{GraphNode, GraphNodeKey, Node};
use crate::make_error;
use eyre::Report;
use itertools::{iproduct, Itertools};
use parking_lot::RwLock;
use std::borrow::BorrowMut;
use std::collections::{HashSet, VecDeque};
use std::fmt::Debug;
use std::io::Write;
use std::sync::Arc;

pub type SafeNode<N> = Arc<RwLock<Node<N>>>;

pub type NodeEdgePair<N, E> = (Arc<RwLock<Node<N>>>, Arc<RwLock<Edge<E>>>);
pub type NodeEdgePayloadPair<N, E> = (Arc<RwLock<N>>, Arc<RwLock<E>>);

/// Represents graph node during forward traversal
#[derive(Debug)]
pub struct GraphNodeForward<'n, N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: &'n mut N,
  pub parents: Vec<NodeEdgePayloadPair<N, E>>,
}

/// Represents graph node during backwards traversal
#[derive(Debug)]
pub struct GraphNodeBackward<'n, N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: &'n mut N,
  pub children: Vec<NodeEdgePayloadPair<N, E>>,
}

/// Represents graph node during safe traversal
#[derive(Debug)]
pub struct GraphNodeSafe<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: Arc<RwLock<N>>,
  pub children: Vec<NodeEdgePair<N, E>>,
  pub parents: Vec<NodeEdgePair<N, E>>,
}

impl<N, E> GraphNodeSafe<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub fn from_node(graph: &Graph<N, E>, node: &Arc<RwLock<Node<N>>>) -> Self {
    let node = node.read();
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();
    let payload = node.payload();
    let parents = graph.parents_of(&node);
    let children = graph.children_of(&node);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      parents,
    }
  }
}

#[derive(Debug)]
pub struct Graph<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  nodes: Vec<Arc<RwLock<Node<N>>>>,
  edges: Vec<Arc<RwLock<Edge<E>>>>,
  roots: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
}

impl<N, E> Graph<N, E>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub const fn new() -> Self {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      roots: vec![],
      leaves: vec![],
    }
  }

  /// Retrieve parent nodes of a given node and the corresponding edges.
  ///
  /// **Returns**: list of pairs `(parent, edge)`, where `parent` is the parent node,
  /// and `edge` is the inbound edge connecting the parent node with the given node.
  pub fn parents_of(&self, node: &Node<N>) -> Vec<NodeEdgePair<N, E>> {
    // Parents are the source nodes of inbound edges
    node
      .inbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .filter_map(|edge| {
        let parent_key = edge.read().source();
        self.get_node(parent_key).map(|parent| (parent, edge))
      })
      .collect_vec()
  }

  /// Retrieve keys of parent nodes of a given node.
  pub fn parent_keys_of(&self, node: &Node<N>) -> Vec<GraphNodeKey> {
    // Parents are the source nodes of inbound edges
    node
      .inbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .map(|edge| edge.read().source())
      .collect_vec()
  }

  /// Retrieve child nodes of a given node and the corresponding edges.
  ///
  /// **Returns**: list of pairs `(child, edge)`, where `child` is the child node,
  /// and `edge` is the outbound edge connecting the given node with the child node.
  pub fn children_of(&self, node: &Node<N>) -> Vec<NodeEdgePair<N, E>> {
    // Children are the target nodes of outbound edges
    node
      .outbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .filter_map(|edge| {
        let child_key = edge.read().target();
        self.get_node(child_key).map(|child| (child, edge))
      })
      .collect_vec()
  }

  /// Retrieve keys of parent nodes of a given node.
  pub fn child_keys_of(&self, node: &Node<N>) -> Vec<GraphNodeKey> {
    // Children are the target nodes of outbound edges
    node
      .outbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .map(|edge| edge.read().target())
      .collect_vec()
  }

  pub fn get_node(&self, index: GraphNodeKey) -> Option<Arc<RwLock<Node<N>>>> {
    self.nodes.get(index.0).map(Arc::clone)
  }

  pub fn get_edge(&self, index: GraphEdgeKey) -> Option<Arc<RwLock<Edge<E>>>> {
    self.edges.get(index.0).map(Arc::clone)
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn for_each<T, F>(&self, f: &mut dyn FnMut(GraphNodeSafe<N, E>)) {
    self
      .nodes
      .iter()
      .for_each(|node| f(GraphNodeSafe::from_node(self, node)));
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E>) -> T,
  {
    self
      .nodes
      .iter()
      .map(|node| f(GraphNodeSafe::from_node(self, node)))
      .collect_vec()
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn filter_map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E>) -> Option<T>,
  {
    self
      .nodes
      .iter()
      .filter_map(|node| f(GraphNodeSafe::from_node(self, node)))
      .collect_vec()
  }

  #[inline]
  pub fn num_nodes(&self) -> usize {
    self.nodes.len()
  }

  #[inline]
  pub fn num_roots(&self) -> usize {
    self.roots.len()
  }

  #[inline]
  pub fn num_leaves(&self) -> usize {
    self.leaves.len()
  }

  #[inline]
  pub fn get_nodes(&self) -> &[Arc<RwLock<Node<N>>>] {
    &self.nodes
  }

  #[inline]
  pub fn get_node_payloads(&self) -> impl Iterator<Item = Arc<RwLock<N>>> + '_ {
    self.nodes.iter().map(|node| node.read().payload())
  }

  #[inline]
  pub fn get_roots(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  #[inline]
  pub fn get_leaves(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self.leaves.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  #[inline]
  pub fn get_edges(&self) -> &[Arc<RwLock<Edge<E>>>] {
    &self.edges
  }

  pub fn add_node(&mut self, node_payload: N) -> GraphNodeKey {
    let node_key = GraphNodeKey(self.nodes.len());
    let node = Arc::new(RwLock::new(Node::new(node_key, node_payload)));
    self.nodes.push(node);
    node_key
  }

  /// Add a new edge to the graph.
  pub fn add_edge(
    &mut self,
    source_key: GraphNodeKey,
    target_key: GraphNodeKey,
    edge_payload: E,
  ) -> Result<(), Report> {
    if source_key == target_key {
      return make_error!(
        "When adding a graph edge {source_key}->{target_key}: Attempted to connect node {source_key} to itself."
      );
    }

    let source_lock = self
      .get_node(source_key)
      .ok_or_else(|| format!("When adding a graph edge {source_key}->{target_key}: Node {source_key} not found."))
      .unwrap();

    let target_lock = self
      .get_node(target_key)
      .ok_or_else(|| format!("When adding a graph edge {source_key}->{target_key}: Node {target_key} not found."))
      .unwrap();

    let edge_key = GraphEdgeKey(self.edges.len());
    let new_edge = Arc::new(RwLock::new(Edge::new(edge_key, source_key, target_key, edge_payload)));

    {
      let (source, target) = (source_lock.read(), target_lock.read());

      let already_connected = source
        .outbound()
        .iter()
        .any(|edge| self.get_edge(*edge).unwrap().read().target() == target.key());

      if already_connected {
        return make_error!("When adding a graph edge {source_key}->{target_key}: Nodes {source_key} and {target_key} are already connected.");
      }

      self.edges.push(Arc::clone(&new_edge));
    }

    {
      let (mut source, mut target) = (source_lock.write(), target_lock.write());
      source.outbound_mut().push(edge_key);
      target.inbound_mut().push(edge_key);
    }

    Ok(())
  }

  pub fn build(&mut self) -> Result<(), Report> {
    self.roots = self
      .nodes
      .iter()
      .filter_map(|node| {
        let node = node.read();
        node.is_root().then(|| node.key())
      })
      .collect_vec();

    self.leaves = self
      .nodes
      .iter()
      .filter_map(|node| {
        let node = node.read();
        node.is_leaf().then(|| node.key())
      })
      .collect_vec();

    Ok(())
  }

  pub fn par_iter_breadth_first_forward<F>(&mut self, explorer: F)
  where
    F: Fn(GraphNodeForward<N, E>) -> GraphTraversalContinuation + Sync + Send,
  {
    let roots = self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    directed_breadth_first_traversal_forward::<N, E, _>(self, roots.as_slice(), |node| {
      let is_leaf = node.is_leaf();
      let is_root = node.is_root();
      let key = node.key();

      let payload = node.payload();
      let mut payload = payload.write();
      let payload = payload.borrow_mut();

      let parents = self
        .parents_of(node)
        .into_iter()
        .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
        .collect_vec();

      explorer(GraphNodeForward {
        is_root,
        is_leaf,
        key,
        payload,
        parents,
      })
    });

    self.reset_nodes();
  }

  pub fn par_iter_breadth_first_backward<F>(&mut self, explorer: F)
  where
    F: Fn(GraphNodeBackward<N, E>) -> GraphTraversalContinuation + Sync + Send,
  {
    let leaves = self.leaves.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    directed_breadth_first_traversal_backward::<N, E, _>(self, leaves.as_slice(), |node| {
      let is_leaf = node.is_leaf();
      let is_root = node.is_root();
      let key = node.key();

      let payload = node.payload();
      let mut payload = payload.write();
      let payload = payload.borrow_mut();

      let children = self
        .children_of(node)
        .into_iter()
        .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
        .collect_vec();

      explorer(GraphNodeBackward {
        is_root,
        is_leaf,
        key,
        payload,
        children,
      })
    });

    self.reset_nodes();
  }

  /// Synchronously traverse graph in depth-first preorder fashion forward (from roots to leaves, along edge directions).
  ///
  /// Guarantees that for each visited node, all of it parents (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_depth_first_preorder_forward(&self, mut explorer: impl FnMut(GraphNodeForward<N, E>)) {
    let mut stack: Vec<GraphNodeKey> = self.roots.clone();
    let mut visited = HashSet::<GraphNodeKey>::new();
    while let Some(key) = stack.pop() {
      if !visited.contains(&key) {
        let node = self.get_node(key).unwrap();
        let node = &node.read();
        stack.extend(self.child_keys_of(node));
        explorer(GraphNodeForward {
          is_root: node.is_root(),
          is_leaf: node.is_leaf(),
          key,
          parents: self
            .parents_of(node)
            .into_iter()
            .map(|(parent, edge)| (parent.read().payload(), edge.read().payload()))
            .collect_vec(),
          payload: &mut node.payload().write(),
        });
        visited.insert(key);
      }
    }
  }

  /// Synchronously traverse graph in breadth-first order forward (from roots to leaves, along edge directions).
  ///
  /// Guarantees that for each visited node, all of it parents (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_breadth_first_forward(&self, mut explorer: impl FnMut(GraphNodeForward<N, E>)) {
    let mut queue: VecDeque<GraphNodeKey> = self.roots.iter().copied().collect();
    let mut visited = HashSet::<GraphNodeKey>::new();
    while let Some(key) = queue.pop_front() {
      if !visited.contains(&key) {
        let node = self.get_node(key).unwrap();
        let node = &node.read();
        queue.extend(self.child_keys_of(node));
        explorer(GraphNodeForward {
          is_root: node.is_root(),
          is_leaf: node.is_leaf(),
          key,
          parents: self
            .parents_of(node)
            .into_iter()
            .map(|(parent, edge)| (parent.read().payload(), edge.read().payload()))
            .collect_vec(),
          payload: &mut node.payload().write(),
        });
        visited.insert(key);
      }
    }
  }

  /// Synchronously traverse graph in breadth-first order backwards (from leaves to roots, against edge directions).
  ///
  /// Guarantees that for each visited node, all of it children (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_breadth_first_reverse(&self, mut explorer: impl FnMut(GraphNodeBackward<N, E>)) {
    let mut queue: VecDeque<GraphNodeKey> = self.leaves.iter().copied().collect();
    let mut visited = HashSet::<GraphNodeKey>::new();
    while let Some(key) = queue.pop_front() {
      if !visited.contains(&key) {
        let node = self.get_node(key).unwrap();
        let node = &node.read();
        queue.extend(self.parent_keys_of(node));
        explorer(GraphNodeBackward {
          is_root: node.is_root(),
          is_leaf: node.is_leaf(),
          key,
          children: self
            .children_of(node)
            .into_iter()
            .map(|(child, edge)| (child.read().payload(), edge.read().payload()))
            .collect_vec(),
          payload: &mut node.payload().write(),
        });
        visited.insert(key);
      }
    }
  }

  /// Returns graph into initial state after traversal
  pub fn reset_nodes(&self) {
    // Mark all nodes as not visited. As a part of traversal all nodes, one by one, marked as visited,
    // to ensure correctness of traversal. Here we reset the "is visited" markers this,
    // to allow for traversals again.
    self.nodes.iter().for_each(|node| node.write().mark_as_not_visited());
  }

  fn print_fake_edges<W: Write>(&self, mut writer: W, nodes: &[&SafeNode<N>]) {
    // Fake edges needed to align a set of nodes beautifully
    let node_keys = nodes.iter().map(|node| node.read().key()).collect_vec();

    let fake_edges = iproduct!(&node_keys, &node_keys)
      .enumerate()
      .map(|(i, (left, right))| {
        let weight = 1000 + i * 100;
        format!("      {left}-> {right} [style=invis, weight={weight}]")
      })
      .join("\n");

    if !fake_edges.is_empty() {
      writeln!(
        writer,
        "\n    // fake edges for alignment of nodes\n    {{\n      rank=same\n{fake_edges}\n    }}"
      )
      .unwrap();
    }
  }

  fn print_node<W: Write>(&self, mut writer: W, node: &SafeNode<N>) {
    writeln!(
      writer,
      "    {} [label = \"({}) {}\"]",
      node.read().key(),
      node.read().key(),
      node.read().payload().read()
    )
    .unwrap();
  }

  /// Print graph nodes.
  fn print_nodes<W: Write>(&self, mut writer: W) {
    let roots = self.nodes.iter().filter(|node| node.read().is_root()).collect_vec();
    let leaves = self.nodes.iter().filter(|node| node.read().is_leaf()).collect_vec();
    let internal = self
      .nodes
      .iter()
      .filter(|node| {
        let node = node.read();
        !node.is_leaf() && !node.is_root()
      })
      .collect_vec();

    writeln!(writer, "\n  subgraph roots {{").unwrap();

    for node in &roots {
      self.print_node(&mut writer, node);
    }

    self.print_fake_edges(&mut writer, &roots);

    writeln!(writer, "  }}\n\n  subgraph internals {{").unwrap();

    for node in &internal {
      self.print_node(&mut writer, node);
    }

    writeln!(writer, "  }}\n\n  subgraph leaves {{").unwrap();

    for node in &leaves {
      self.print_node(&mut writer, node);
    }

    self.print_fake_edges(&mut writer, &leaves);

    writeln!(writer, "  }}").unwrap();
  }

  /// Print graph edges.
  fn print_edges<W: Write>(&self, mut writer: W) {
    self.nodes.iter().for_each(&mut |node: &Arc<RwLock<Node<N>>>| {
      for edge_key in node.read().outbound() {
        let edge = self.get_edge(*edge_key).unwrap();
        let payload = edge.read().payload();
        let payload = payload.read();
        writeln!(
          writer,
          "  {} -> {} [xlabel = \"{}\", weight=\"{}\"]",
          edge.read().source(),
          edge.read().target(),
          payload,
          payload.weight()
        )
        .unwrap();
      }
    });
  }

  /// Print graph in .dot format.
  pub fn print_graph<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
    write!(
      writer,
      r#"
digraph Phylogeny {{
  graph [rankdir=LR, overlap=scale, splines=ortho, nodesep=1.0, ordering=out];
  edge  [overlap=scale];
  node  [shape=box];
"#
    )?;
    self.print_nodes(&mut writer);
    writeln!(writer)?;
    self.print_edges(&mut writer);
    writeln!(writer, "}}")?;
    Ok(())
  }
}
