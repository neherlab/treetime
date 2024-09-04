use crate::graph::breadth_first::{
  directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward, GraphTraversalContinuation,
};
use crate::graph::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::node::{GraphNode, GraphNodeKey, Node};
use crate::utils::container::get_exactly_one_mut;
use crate::{make_error, make_internal_error, make_internal_report, make_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::lock_api::{ArcRwLockReadGuard, ArcRwLockWriteGuard};
use parking_lot::{RawRwLock, RwLock};
use serde::{Deserialize, Serialize};
use std::collections::{HashSet, VecDeque};
use std::fmt::Debug;
use std::sync::Arc;
use traversal::{Bft, DftPre};

pub type SafeNode<N> = Arc<RwLock<Node<N>>>;
pub type SafeNodeRef<N> = ArcRwLockReadGuard<RawRwLock, Node<N>>;
pub type SafeNodeRefMut<N> = ArcRwLockWriteGuard<RawRwLock, Node<N>>;

pub type SafeEdge<E> = Arc<RwLock<Edge<E>>>;
pub type SafeEdgeRef<E> = ArcRwLockReadGuard<RawRwLock, Edge<E>>;
pub type SafeEdgeRefMut<E> = ArcRwLockWriteGuard<RawRwLock, Edge<E>>;

pub type SafeNodePayloadRef<N> = ArcRwLockReadGuard<RawRwLock, N>;
pub type SafeNodePayloadRefMut<N> = ArcRwLockWriteGuard<RawRwLock, N>;

pub type SafeEdgePayloadRef<E> = ArcRwLockReadGuard<RawRwLock, E>;
pub type SafeEdgePayloadRefMut<E> = ArcRwLockWriteGuard<RawRwLock, E>;

pub type NodeEdgePair<N, E> = (Arc<RwLock<Node<N>>>, Arc<RwLock<Edge<E>>>);
pub type NodeEdgePayloadPair<N, E> = (Arc<RwLock<N>>, Arc<RwLock<E>>);
pub type NodeRefEdgeRefMutPayloadPair<N, E> = (SafeNodePayloadRef<N>, SafeEdgePayloadRefMut<E>);

/// Represents graph node during forward traversal
#[must_use]
#[derive(Debug)]
pub struct GraphNodeForward<'a, N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  graph: &'a Graph<N, E, D>,
  node: &'a Node<N>,
}

impl<'g, N, E, D> GraphNodeForward<'g, N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn new(graph: &'g Graph<N, E, D>, node: &'g Node<N>) -> Self {
    Self { graph, node }
  }

  #[must_use]
  pub fn is_root(&self) -> bool {
    self.node.is_root()
  }

  #[must_use]
  pub fn is_leaf(&self) -> bool {
    self.node.is_leaf()
  }

  #[must_use]
  pub fn key(&self) -> GraphNodeKey {
    self.node.key()
  }

  /// Mutable reference to node's payload
  pub fn payload(&self) -> SafeNodePayloadRef<N> {
    self.node.payload().read_arc()
  }

  /// Mutable reference to node's payload
  pub fn payload_mut(&mut self) -> SafeNodePayloadRefMut<N> {
    self.node.payload().write_arc()
  }

  /// Iterator of pairs (parent, edge), where:
  ///  * parent - reference to parent node payload
  ///  * edge - mutable reference to payload of the edge connecting this node and the parent node
  #[must_use]
  pub fn parents(&self) -> impl DoubleEndedIterator<Item = NodeRefEdgeRefMutPayloadPair<N, E>> + '_ {
    self.graph.parents_of(self.node).map(|(parent, edge)| {
      (
        parent.read_arc().payload().read_arc(),
        edge.read_arc().payload().write_arc(),
      )
    })
  }

  /// Same as .parents(), but:
  ///  * if 1 parent: returns only first parent
  ///  * if no parents: returns None
  ///  * if more than 2 parent: returns error
  pub fn get_at_most_one_parent(&self) -> Result<Option<NodeRefEdgeRefMutPayloadPair<N, E>>, Report> {
    self.parents().at_most_one().map_err(|i| {
      make_report!(
        "Multiple parent nodes are not supported yet, but node {} has {} parents",
        self.node.key(),
        i.count()
      )
    })
  }

  /// Same as .parents(), but:
  ///  * if 1 parent: returns only first parent
  ///  * otherwise returns error
  pub fn get_exactly_one_parent(&self) -> Result<NodeRefEdgeRefMutPayloadPair<N, E>, Report> {
    self.parents().exactly_one().map_err(|i| {
      make_report!(
        "Multiple parent nodes are not supported yet, but node {} has {} parents",
        self.node.key(),
        i.count()
      )
    })
  }

  /// Iterator of mutable references to payloads of edges connecting this node and its children
  pub fn child_edges(&self) -> impl DoubleEndedIterator<Item = SafeEdgePayloadRefMut<E>> + 'g {
    self
      .graph
      .children_of(self.node)
      .map(|(_, edge)| edge.write_arc().payload().write_arc())
  }
}

/// Represents graph node during backwards traversal
#[must_use]
#[derive(Debug)]
pub struct GraphNodeBackward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: SafeNodePayloadRefMut<N>,
  pub children: Vec<NodeEdgePayloadPair<N, E>>,
  pub parent_edges: Vec<SafeEdgePayloadRefMut<E>>,
  pub data: Arc<RwLock<D>>,
}

impl<N, E, D> GraphNodeBackward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn new(graph: &Graph<N, E, D>, node: &Node<N>) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let payload = node.payload().write_arc();

    let children = graph
      .children_of(node)
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    let parent_edges = graph
      .parents_of(node)
      .map(|(_, edge)| edge.write_arc().payload().write_arc())
      .collect_vec();

    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      parent_edges,
      data,
    }
  }

  pub fn get_exactly_one_parent_edge(&mut self) -> Result<&mut SafeEdgePayloadRefMut<E>, Report> {
    get_exactly_one_mut(&mut self.parent_edges).wrap_err("Nodes with multiple parents are not yet supported")
  }
}

/// Represents graph node during safe traversal
#[derive(Debug)]
pub struct GraphNodeSafe<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: Arc<RwLock<N>>,
  pub children: Vec<NodeEdgePair<N, E>>,
  pub parents: Vec<NodeEdgePair<N, E>>,
  pub data: Arc<RwLock<D>>,
}

impl<N, E, D> GraphNodeSafe<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn from_node(graph: &Graph<N, E, D>, node: &Arc<RwLock<Node<N>>>) -> Self {
    let node = node.read();
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();
    let payload = node.payload();
    let parents = graph.parents_of(&node).collect_vec();
    let children = graph.children_of(&node).collect_vec();
    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      parents,
      data,
    }
  }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<N, E, D = ()>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  nodes: Vec<Option<Arc<RwLock<Node<N>>>>>,
  edges: Vec<Option<Arc<RwLock<Edge<E>>>>>,
  roots: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
  data: Arc<RwLock<D>>,
}

impl<N, E, D> Graph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn new() -> Self
  where
    D: Default,
  {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      roots: vec![],
      leaves: vec![],
      data: Arc::new(RwLock::new(D::default())),
    }
  }

  pub fn with_data(data: D) -> Self {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      roots: vec![],
      leaves: vec![],
      data: Arc::new(RwLock::new(data)),
    }
  }

  pub fn data(&self) -> Arc<RwLock<D>> {
    Arc::clone(&self.data)
  }

  pub fn set_data(&mut self, data: D) {
    self.data = Arc::new(RwLock::new(data));
  }

  /// Retrieve parent nodes of a given node and the corresponding edges.
  ///
  /// **Returns**: list of pairs `(parent, edge)`, where `parent` is the parent node,
  /// and `edge` is the inbound edge connecting the parent node with the given node.
  pub fn parents_of<'a>(&'a self, node: &'a Node<N>) -> impl DoubleEndedIterator<Item = NodeEdgePair<N, E>> + 'a {
    // Parents are the source nodes of inbound edges
    node
      .inbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .filter_map(|edge| {
        let parent_key = edge.read_arc().source();
        self.get_node(parent_key).map(|parent| (parent, edge))
      })
  }

  /// Retrieve keys of parent nodes of a given node.
  pub fn parent_keys_of<'a>(&'a self, node: &'a Node<N>) -> impl DoubleEndedIterator<Item = GraphNodeKey> + 'a {
    // Parents are the source nodes of inbound edges
    node
      .inbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .map(|edge| edge.read_arc().source())
  }

  pub fn exactly_one_parent_of(&self, node: &Node<N>) -> Result<NodeEdgePair<N, E>, Report> {
    self.one_parent_of(node)?.ok_or_else(|| {
      make_internal_report!(
        "No parents found for node {} (context: is_root={} is_leaf={})",
        node.key(),
        node.is_root(),
        node.is_leaf()
      )
    })
  }

  pub fn one_parent_of(&self, node: &Node<N>) -> Result<Option<NodeEdgePair<N, E>>, Report> {
    let parents = self.parents_of(node).collect_vec();

    if parents.is_empty() {
      return Ok(None);
    }

    if parents.len() > 1 {
      return make_internal_error!(
        "Only trees with exactly one parent per node are currently supported, but node '{}' has {} parents",
        node.key(),
        self.roots.len()
      );
    }

    Ok(Some(parents[0].clone()))
  }

  /// Retrieve child nodes of a given node and the corresponding edges.
  ///
  /// **Returns**: list of pairs `(child, edge)`, where `child` is the child node,
  /// and `edge` is the outbound edge connecting the given node with the child node.
  pub fn children_of<'a>(&'a self, node: &'a Node<N>) -> impl DoubleEndedIterator<Item = NodeEdgePair<N, E>> + 'a {
    node
      .outbound()
      .iter()
      .filter_map(|edge_key| self.get_edge(*edge_key))
      .map(|edge| {
        let child_key = edge.read_arc().target();
        let child = self.get_node(child_key).expect("Edge not found");
        (child, edge)
      })
  }

  pub fn children_of_by_key<'a>(
    &'a self,
    node_key: GraphNodeKey,
  ) -> impl DoubleEndedIterator<Item = NodeEdgePair<N, E>> + 'a {
    // Children are the target nodes of outbound edges
    let outbound = self
      .get_node(node_key)
      .expect("Node not found")
      .read_arc()
      .outbound()
      .to_owned();

    outbound
      .into_iter()
      .filter_map(|edge_key| self.get_edge(edge_key))
      .filter_map(|edge| {
        let child_key = edge.read().target();
        self.get_node(child_key).map(|child| (child, edge))
      })
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
    self.nodes.get(index.as_usize())?.as_ref().map(Arc::clone)
  }

  pub fn get_edge(&self, index: GraphEdgeKey) -> Option<Arc<RwLock<Edge<E>>>> {
    self.edges.get(index.as_usize())?.as_ref().map(Arc::clone)
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn for_each<T, F>(&self, f: &mut dyn FnMut(GraphNodeSafe<N, E, D>)) {
    self
      .nodes
      .iter()
      .filter_map(Option::as_ref)
      .for_each(|node| f(GraphNodeSafe::from_node(self, node)));
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E, D>) -> T,
  {
    self
      .nodes
      .iter()
      .filter_map(Option::as_ref)
      .map(|node| f(GraphNodeSafe::from_node(self, node)))
      .collect_vec()
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn filter_map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E, D>) -> Option<T>,
  {
    self
      .nodes
      .iter()
      .filter_map(Option::as_ref)
      .filter_map(|node| f(GraphNodeSafe::from_node(self, node)))
      .collect_vec()
  }

  pub fn num_nodes(&self) -> usize {
    self.nodes.len()
  }

  pub fn num_roots(&self) -> usize {
    self.roots.len()
  }

  pub fn num_leaves(&self) -> usize {
    self.leaves.len()
  }

  // All nodes
  pub fn get_nodes(&self) -> Vec<SafeNode<N>> {
    self.nodes.iter().filter_map(Option::as_ref).cloned().collect()
  }

  pub fn get_node_payloads(&self) -> impl Iterator<Item = Arc<RwLock<N>>> + '_ {
    self
      .nodes
      .iter()
      .filter_map(|node| node.as_ref().map(|n| Arc::clone(&n.read_arc().payload())))
  }

  pub fn get_exactly_one_root(&self) -> Result<Arc<RwLock<Node<N>>>, Report> {
    let roots = self.get_roots();
    if roots.len() != 1 {
      make_internal_error!(
        "Only trees with exactly one root are currently supported, but found {} roots",
        self.roots.len()
      )
    } else {
      Ok(Arc::clone(&roots[0]))
    }
  }

  // All nodes having no parents
  pub fn get_roots(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  pub fn get_root_payloads(&self) -> impl Iterator<Item = Arc<RwLock<N>>> + '_ {
    self
      .roots
      .iter()
      .filter_map(|idx| self.get_node(*idx))
      .map(|node| node.read().payload())
  }

  // All nodes having no children
  pub fn get_leaves(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self.leaves.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  // All nodes which are not leaves
  pub fn get_internal_nodes(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self
      .get_nodes()
      .iter()
      .filter(|node| !node.read_arc().is_leaf())
      .map(Arc::clone)
      .collect_vec()
  }

  // All nodes which are neither leaves nor roots (i.e. internal which are not roots)
  pub fn get_inner_nodes(&self) -> Vec<Arc<RwLock<Node<N>>>> {
    self
      .get_nodes()
      .iter()
      .filter(|node| !node.read_arc().is_leaf() && !node.read_arc().is_root())
      .map(Arc::clone)
      .collect_vec()
  }

  pub fn get_edges(&self) -> Vec<Arc<RwLock<Edge<E>>>> {
    self.edges.iter().filter_map(Option::as_ref).cloned().collect()
  }

  pub fn add_node(&mut self, node_payload: N) -> GraphNodeKey {
    let node_key = GraphNodeKey(self.nodes.len());
    let node = Arc::new(RwLock::new(Node::new(node_key, node_payload)));
    self.nodes.push(Some(node));
    node_key
  }

  pub fn remove_node(&mut self, node_key: GraphNodeKey) -> Result<(), Report> {
    // Remove edges which refer to this node, if any
    self
      .edges
      .iter_mut()
      // Find edges which refer to this node
      .filter(|e| {
        e.as_ref().map_or(false, |e| {
          let e = e.read();
          e.source() == node_key || e.target() == node_key
        })
      })
      // Remove these edges
      .for_each(|edge| *edge = None);

    // Remove the node itself, if present
    if let Some(node) = self.nodes.get_mut(node_key.as_usize()) {
      *node = None;
    }

    Ok(())
  }

  /// Add a new edge to the graph.
  pub fn add_edge(
    &mut self,
    source_key: GraphNodeKey,
    target_key: GraphNodeKey,
    edge_payload: E,
  ) -> Result<GraphEdgeKey, Report> {
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

      self.edges.push(Some(Arc::clone(&new_edge)));
    }

    {
      let (mut source, mut target) = (source_lock.write(), target_lock.write());
      source.outbound_mut().push(edge_key);
      target.inbound_mut().push(edge_key);
    }

    Ok(edge_key)
  }

  pub fn remove_edge(&mut self, edge_key: GraphEdgeKey) -> Result<(), Report> {
    // Remove the edge key from inbound/outbound lists of nodes
    self.nodes.iter_mut().for_each(|node| {
      if let Some(node) = node {
        let mut node_locked = node.write_arc();
        node_locked.outbound_mut().retain(|&e| e != edge_key);
        node_locked.inbound_mut().retain(|&e| e != edge_key);
      }
    });

    // Remove the edge itself
    if let Some(edge) = self.edges.get_mut(edge_key.as_usize()) {
      *edge = None;
    }

    Ok(())
  }

  pub fn build(&mut self) -> Result<(), Report> {
    self.roots = self
      .nodes
      .iter()
      .filter_map(|node_option| {
        node_option
          .as_ref()
          .and_then(|node| node.read().is_root().then(|| node.read().key()))
      })
      .collect();

    self.leaves = self
      .nodes
      .iter()
      .filter_map(|node_option| {
        node_option
          .as_ref()
          .and_then(|node| node.read().is_leaf().then(|| node.read().key()))
      })
      .collect();

    Ok(())
  }

  pub fn par_iter_breadth_first_forward<F>(&self, explorer: F)
  where
    F: Fn(GraphNodeForward<N, E, D>) -> GraphTraversalContinuation + Sync + Send,
  {
    let roots = self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    directed_breadth_first_traversal_forward::<N, E, D, _>(self, roots.as_slice(), |node| {
      explorer(GraphNodeForward::new(self, node))
    });
  }

  pub fn par_iter_breadth_first_backward<F>(&self, explorer: F)
  where
    F: Fn(GraphNodeBackward<N, E, D>) -> GraphTraversalContinuation + Sync + Send,
  {
    let leaves = self
      .leaves
      .iter()
      .filter_map(|idx| self.get_node(*idx))
      .rev()
      .collect_vec();

    directed_breadth_first_traversal_backward::<N, E, D, _>(self, leaves.as_slice(), |node| {
      explorer(GraphNodeBackward::new(self, node))
    });
  }

  /// Synchronously traverse graph in depth-first preorder fashion forward (from roots to leaves, along edge directions).
  ///
  /// Guarantees that for each visited node, all of it parents (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_depth_first_preorder_forward(&self, mut explorer: impl FnMut(GraphNodeForward<N, E, D>)) {
    let root = self.get_exactly_one_root().unwrap();
    let mut stack = Vec::from([(Arc::clone(&root), None)]);
    while let Some((current_node, _current_edge)) = stack.pop() {
      let current_node = current_node.read_arc();
      explorer(GraphNodeForward::new(self, &current_node));
      for (child, edge) in self.children_of(&current_node).rev() {
        stack.push((child, Some(edge)));
      }
    }
  }

  pub fn iter_depth_first_preorder_forward_2(&self, mut explorer: impl FnMut(&SafeNode<N>)) {
    let root = self.get_exactly_one_root().unwrap();
    DftPre::new(&root, |node| self.iter_children_arc(node)).for_each(move |(_, node)| {
      explorer(node);
    });
  }

  /// Synchronously traverse graph in depth-first postorder fashion forward (from roots to leaves, along edge directions).
  ///
  /// Guarantees that for each visited node, all of it children (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_depth_first_postorder_forward(&self, mut explorer: impl FnMut(GraphNodeBackward<N, E, D>)) {
    let root = self.get_exactly_one_root().unwrap();
    let mut stack = Vec::new();
    let mut visited = HashSet::new();
    stack.push((Arc::clone(&root), None));
    while let Some((current_node, _)) = stack.pop() {
      let node_key = current_node.read_arc().key();
      if !visited.contains(&node_key) {
        visited.insert(node_key);
        stack.push((Arc::clone(&current_node), None));
        let current_node = current_node.read_arc();
        let children = self.children_of(&current_node).rev();
        for (child, edge) in children {
          let child_key = child.read_arc().key();
          if !visited.contains(&child_key) {
            stack.push((child, Some(edge)));
          }
        }
      } else {
        explorer(GraphNodeBackward::new(self, &current_node.read_arc()));
      }
    }
  }

  /// Synchronously traverse graph in breadth-first order forward (from roots to leaves, along edge directions).
  ///
  /// Guarantees that for each visited node, all of it parents (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_breadth_first_forward(&self, mut explorer: impl FnMut(GraphNodeForward<N, E, D>)) {
    let root = self.get_exactly_one_root().unwrap();
    let mut queue = VecDeque::new();
    queue.push_back(Arc::clone(&root));

    while let Some(current_node) = queue.pop_front() {
      explorer(GraphNodeForward::new(self, &current_node.read_arc()));
      for (child, _) in self.children_of(&current_node.read_arc()) {
        queue.push_back(child);
      }
    }
  }

  /// Synchronously traverse graph in breadth-first order backwards (from leaves to roots, against edge directions).
  ///
  /// Guarantees that for each visited node, all of it children (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_breadth_first_reverse(&self, mut explorer: impl FnMut(GraphNodeBackward<N, E, D>)) {
    let root = self.get_exactly_one_root().unwrap();
    Bft::new(&root, |node| self.iter_children_arc(node))
      .collect_vec() // HACK: how to reverse without collecting?
      .into_iter()
      .rev()
      .for_each(move |(_, node)| {
        explorer(GraphNodeBackward::new(self, &node.write()));
      });
  }

  fn iter_children_arc(&self, node: &Arc<RwLock<Node<N>>>) -> impl Iterator<Item = &Arc<RwLock<Node<N>>>> {
    let child_keys = self.child_keys_of(&*node.read());
    self.nodes.iter().filter_map(move |node| {
      node
        .as_ref()
        .and_then(|node| child_keys.contains(&node.read_arc().key()).then_some(node))
    })
  }

  /// Find nodes on the path from root to a given node
  pub fn path_from_root_to_node(&self, node_key: GraphNodeKey) -> Result<Vec<SafeNode<N>>, Report> {
    let mut node = self
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;

    let mut path = vec![Arc::clone(&node)];
    loop {
      match self.one_parent_of(&node.read_arc())? {
        None => break,
        Some((parent, _)) => {
          path.push(Arc::clone(&parent));
          node = parent;
        }
      }
    }

    let path = path.into_iter().rev().collect_vec();
    Ok(path)
  }

  #[allow(clippy::type_complexity)]
  pub fn path_from_node_to_node(
    &self,
    start: GraphNodeKey,
    finish: GraphNodeKey,
  ) -> Result<Vec<(SafeNode<N>, Option<SafeEdge<E>>)>, Report> {
    let mut node = self
      .get_node(start)
      .ok_or_else(|| make_internal_report!("Node not found on the graph: {start}"))?;

    let mut path = vec![(Arc::clone(&node), None)];
    loop {
      let parent = self.one_parent_of(&node.read_arc())?;

      match parent {
        None => {
          return make_internal_error!("When searching path from starting node {start} to destination node {finish}: reached root node without finding the destination")
        }

        Some((parent, edge)) => {
          path.push((Arc::clone(&node), Some(Arc::clone(&edge))));

          if parent.read_arc().key() == finish {
            break;
          }

          node = parent;
        }
      }
    }

    Ok(path)
  }

  pub fn is_root(&self, key: GraphNodeKey) -> bool {
    self.roots.contains(&key)
  }

  pub fn is_leaf(&self, key: GraphNodeKey) -> bool {
    self.leaves.contains(&key)
  }
}

#[cfg(test)]
pub mod tests {
  use std::collections::BTreeMap;

  use pretty_assertions::assert_eq;

  use crate::graph::edge::Weighted;
  use crate::graph::node::Named;
  use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
  use crate::io::nwk::{format_weight, nwk_read_str, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};

  use super::*;

  pub type TestGraph = Graph<TestEdge, TestNode, ()>;

  #[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
  pub struct TestNode(pub Option<String>);

  impl GraphNode for TestNode {}

  impl Named for TestNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      self.0.as_deref()
    }
    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.0 = name.map(|n| n.as_ref().to_owned());
    }
  }

  impl NodeFromNwk for TestNode {
    fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
      Ok(Self(name.map(|n| n.as_ref().to_owned())))
    }
  }

  impl NodeToNwk for TestNode {
    fn nwk_name(&self) -> Option<impl AsRef<str>> {
      self.0.as_deref()
    }
  }

  impl NodeToGraphviz for TestNode {
    fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
      self.0.as_deref()
    }
  }

  #[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
  pub struct TestEdge(pub Option<f64>);

  impl GraphEdge for TestEdge {}

  impl Weighted for TestEdge {
    fn weight(&self) -> Option<f64> {
      self.0
    }
    fn set_weight(&mut self, weight: Option<f64>) {
      self.0 = weight;
    }
  }

  impl EdgeFromNwk for TestEdge {
    fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
      Ok(Self(weight))
    }
  }

  impl EdgeToNwk for TestEdge {
    fn nwk_weight(&self) -> Option<f64> {
      self.0
    }
  }

  impl EdgeToGraphViz for TestEdge {
    fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
      self.0.map(|weight| format_weight(weight, &NwkWriteOptions::default()))
    }

    fn to_graphviz_weight(&self) -> Option<f64> {
      self.0
    }
  }

  #[test]
  fn test_traversal_serial_depth_first_preorder_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_depth_first_preorder_forward(|node| {
      actual.push(node.payload().name().unwrap().as_ref().to_owned());
    });

    assert_eq!(vec!["root", "AB", "A", "B", "CD", "C", "D"], actual);

    Ok(())
  }

  #[test]
  fn test_traversal_serial_depth_first_postorder_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_depth_first_postorder_forward(|node| {
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
    });

    assert_eq!(vec!["A", "B", "AB", "C", "D", "CD", "root"], actual);

    Ok(())
  }

  #[test]
  fn test_traversal_serial_breadth_first_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_breadth_first_forward(|node| {
      actual.push(node.payload().name().unwrap().as_ref().to_owned());
    });

    assert_eq!(vec!["root", "AB", "CD", "A", "B", "C", "D"], actual);

    Ok(())
  }

  #[test]
  fn test_traversal_serial_breadth_first_reverse() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_breadth_first_reverse(|node| {
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
    });

    assert_eq!(vec!["D", "C", "B", "A", "CD", "AB", "root"], actual);

    Ok(())
  }

  #[test]
  fn test_traversal_parallel_breadth_first_forward() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_forward(|node| {
      actual
        .write_arc()
        .push(node.payload().name().unwrap().as_ref().to_owned());
      GraphTraversalContinuation::Continue
    });

    assert_eq!(&vec!["root", "AB", "CD", "A", "B", "C", "D"], &*actual.read());

    Ok(())
  }

  #[test]
  fn test_traversal_parallel_breadth_first_backward() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_backward(|node| {
      actual
        .write_arc()
        .push(node.payload.name().unwrap().as_ref().to_owned());
      GraphTraversalContinuation::Continue
    });

    assert_eq!(&vec!["D", "C", "B", "A", "CD", "AB", "root"], &*actual.read());

    Ok(())
  }
}
