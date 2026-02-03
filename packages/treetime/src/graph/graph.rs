use crate::graph::breadth_first::{
  GraphTraversalContinuation, directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward,
};
use crate::graph::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph::node::{GraphNode, GraphNodeKey, Node};
use crate::{make_error, make_internal_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::lock_api::{ArcRwLockReadGuard, ArcRwLockWriteGuard};
use parking_lot::{RawRwLock, RwLock};
use serde::{Deserialize, Serialize};
use std::collections::{HashSet, VecDeque};
use std::fmt::Debug;
use std::sync::Arc;
use traversal::{Bft, DftPre};
use treetime_utils::container::{get_exactly_one, get_exactly_one_mut};
use treetime_utils::mutex::unwrap_arc_rwlock;

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

/// Represents graph node during forward traversal
#[must_use]
#[derive(Debug)]
pub struct GraphNodeForward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: GraphNodeKey,
  pub payload: SafeNodePayloadRefMut<N>,
  pub parents: Vec<NodeEdgePayloadPair<N, E>>,
  pub parent_keys: Vec<(GraphNodeKey, GraphEdgeKey)>,
  pub child_edges: Vec<SafeEdgePayloadRefMut<E>>,
  pub child_edge_keys: Vec<GraphEdgeKey>,
  pub data: Arc<RwLock<D>>,
}

impl<N, E, D> GraphNodeForward<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub fn new(graph: &Graph<N, E, D>, node: &Node<N>) -> Self {
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

    let parent_keys = graph
      .parents_of(node)
      .iter()
      .map(|(node, edge)| (node.read_arc().key(), edge.read_arc().key()))
      .collect_vec();

    let child_edge_keys = graph
      .children_of(node)
      .iter()
      .map(|(_, edge)| edge.read_arc().key())
      .collect_vec();

    let payload = node.payload().write_arc();

    let parents = graph
      .parents_of(node)
      .iter()
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    let child_edges = graph
      .children_of(node)
      .iter()
      .map(|(_, edge)| edge.write_arc().payload().write_arc())
      .collect_vec();

    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      parents,
      parent_keys,
      child_edges,
      child_edge_keys,
      data,
    }
  }

  pub fn get_exactly_one_parent(&self) -> Result<NodeEdgePayloadPair<N, E>, Report> {
    get_exactly_one(&self.parents)
      .cloned()
      .wrap_err("Nodes with multiple parents are not yet supported")
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
  pub child_keys: Vec<(GraphNodeKey, GraphEdgeKey)>,
  pub parent_edges: Vec<SafeEdgePayloadRefMut<E>>,
  pub parent_edge_keys: Vec<GraphEdgeKey>,
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

    let child_keys = graph
      .children_of(node)
      .iter()
      .map(|(node, edge)| (node.read_arc().key(), edge.read_arc().key()))
      .collect_vec();

    let parent_edge_keys = graph
      .parents_of(node)
      .iter()
      .map(|(_, edge)| edge.read_arc().key())
      .collect_vec();

    let payload = node.payload().write_arc();

    let children = graph
      .children_of(node)
      .iter()
      .map(|(node, edge)| (node.read().payload(), edge.read().payload()))
      .collect_vec();

    let parent_edges = graph
      .parents_of(node)
      .iter()
      .map(|(_, edge)| edge.write_arc().payload().write_arc())
      .collect_vec();

    let data = Arc::clone(&graph.data);

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      child_keys,
      parent_edges,
      parent_edge_keys,
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
    let parents = graph.parents_of(&node);
    let children = graph.children_of(&node);
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
    let parents = self.parents_of(node).into_iter().collect_vec();

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
    self.nodes.get(index.as_usize())?.as_ref().map(Arc::clone)
  }

  pub fn find_node<F>(&self, predicate: F) -> Option<GraphNodeKey>
  where
    F: Fn(&N) -> bool,
  {
    self.nodes.iter().find_map(|node| {
      node.as_ref().and_then(|node| {
        let node_guard = node.read_arc();
        let payload = node_guard.payload().read_arc();
        predicate(&payload).then_some(node_guard.key())
      })
    })
  }

  pub fn get_edge(&self, index: GraphEdgeKey) -> Option<Arc<RwLock<Edge<E>>>> {
    self.edges.get(index.as_usize())?.as_ref().map(Arc::clone)
  }

  pub fn get_source_node_key(&self, edge_key: GraphEdgeKey) -> Result<GraphNodeKey, Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    Ok(edge.read_arc().source())
  }

  pub fn get_target_node_key(&self, edge_key: GraphEdgeKey) -> Result<GraphNodeKey, Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    Ok(edge.read_arc().target())
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

  #[allow(clippy::needless_collect)]
  pub fn remove_node(&mut self, node_key: GraphNodeKey) -> Result<(Node<N>, Vec<Edge<E>>), Report> {
    let edges_to_remove: Vec<GraphEdgeKey> = self
      .edges
      .iter()
      .filter_map(|e| {
        e.as_ref().and_then(|e| {
          let e = e.read();
          (e.source() == node_key || e.target() == node_key).then_some(e.key())
        })
      })
      .collect();

    let removed_edges = edges_to_remove
      .into_iter()
      .map(|edge_key| -> Result<Edge<E>, Report> { self.remove_edge(edge_key) })
      .collect::<Result<Vec<_>, _>>()?;

    let removed_node = self
      .nodes
      .get_mut(node_key.as_usize())
      .and_then(|node| node.take().map(unwrap_arc_rwlock))
      .transpose()
      .wrap_err_with(|| format!("When removing node: {node_key}"))?
      .ok_or_else(|| make_internal_report!("Attempted to remove non-existent node: {node_key}"))?;

    Ok((removed_node, removed_edges))
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
        return make_error!(
          "When adding a graph edge {source_key}->{target_key}: Nodes {source_key} and {target_key} are already connected."
        );
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

  pub fn remove_edge(&mut self, edge_key: GraphEdgeKey) -> Result<Edge<E>, Report> {
    // Remove the edge key from inbound/outbound lists of nodes
    self.nodes.iter_mut().for_each(|node| {
      if let Some(node) = node {
        let mut node_locked = node.write_arc();
        node_locked.outbound_mut().retain(|&e| e != edge_key);
        node_locked.inbound_mut().retain(|&e| e != edge_key);
      }
    });

    // Remove the edge itself
    self
      .edges
      .get_mut(edge_key.as_usize())
      .and_then(|edge_slot| edge_slot.take().map(unwrap_arc_rwlock))
      .transpose()
      .wrap_err_with(|| format!("When removing edge: {edge_key}"))?
      .ok_or_else(|| make_internal_report!("Attempted to remove non-existent edge: {edge_key}"))
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

    self.reset_nodes();
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

    self.reset_nodes();
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
      for (child, edge) in self.children_of(&current_node).into_iter().rev() {
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
      if visited.insert(node_key) {
        stack.push((Arc::clone(&current_node), None));
        let children = self.children_of(&current_node.read_arc()).into_iter().rev();
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
      let children = self.children_of(&current_node.read_arc());
      for (child, _) in children {
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
        },
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
          return make_internal_error!(
            "When searching path from starting node {start} to destination node {finish}: reached root node without finding the destination"
          );
        },

        Some((parent, edge)) => {
          path.push((Arc::clone(&node), Some(Arc::clone(&edge))));

          if parent.read_arc().key() == finish {
            break;
          }

          node = parent;
        },
      }
    }

    Ok(path)
  }

  /// Returns graph into initial state after traversal
  pub fn reset_nodes(&self) {
    // Mark all nodes as not visited. As a part of traversal all nodes, one by one, marked as visited,
    // to ensure correctness of traversal. Here we reset the "is visited" markers this,
    // to allow for traversals again.
    self
      .nodes
      .iter()
      .filter_map(|node| node.as_ref())
      .for_each(|node| node.write().mark_as_not_visited());
  }

  pub fn is_root(&self, key: GraphNodeKey) -> bool {
    self.roots.contains(&key)
  }

  pub fn is_leaf(&self, node_key: GraphNodeKey) -> bool {
    self.leaves.contains(&node_key)
  }

  pub fn is_internal(&self, node_key: GraphNodeKey) -> bool {
    !self.is_leaf(node_key) && !self.is_root(node_key)
  }

  /// Number of outbound edges (children) of a node
  pub fn degree_out(&self, key: GraphNodeKey) -> Result<usize, Report> {
    self
      .get_node(key)
      .map(|node| node.read_arc().degree_out())
      .ok_or_else(|| make_internal_report!("Node not found: {key}"))
  }

  /// Number of inbound edges (parents) of a node
  pub fn degree_in(&self, key: GraphNodeKey) -> Result<usize, Report> {
    self
      .get_node(key)
      .map(|node| node.read_arc().degree_in())
      .ok_or_else(|| make_internal_report!("Node not found: {key}"))
  }

  pub fn has_parents(&self, node_key: GraphNodeKey) -> bool {
    self
      .get_node(node_key)
      .is_some_and(|node| node.read_arc().has_parents())
  }

  pub fn has_one_parent(&self, node_key: GraphNodeKey) -> bool {
    self
      .get_node(node_key)
      .is_some_and(|node| node.read_arc().has_one_parent())
  }

  pub fn has_at_most_one_parent(&self, node_key: GraphNodeKey) -> bool {
    self
      .get_node(node_key)
      .is_some_and(|node| node.read_arc().has_at_most_one_parent())
  }

  pub fn has_children(&self, node_key: GraphNodeKey) -> bool {
    self
      .get_node(node_key)
      .is_some_and(|node| node.read_arc().has_children())
  }

  pub fn has_one_child(&self, node_key: GraphNodeKey) -> bool {
    self
      .get_node(node_key)
      .is_some_and(|node| node.read_arc().has_one_child())
  }

  pub fn has_at_most_one_child(&self, node_key: GraphNodeKey) -> bool {
    self
      .get_node(node_key)
      .is_some_and(|node| node.read_arc().has_at_most_one_child())
  }

  /// Returns the inbound edge from this node's parent, if any (None for roots).
  pub fn parent_inbound_edge(&self, key: GraphNodeKey) -> Result<Option<GraphEdgeKey>, Report> {
    let node = self
      .get_node(key)
      .ok_or_else(|| make_internal_report!("Node not found: {key}"))?;
    Ok(node.read_arc().inbound().first().copied())
  }

  #[allow(clippy::type_complexity)]
  pub fn collapse_edge(&mut self, edge_key: GraphEdgeKey) -> Result<(Node<N>, Edge<E>, Vec<SafeEdge<E>>), Report>
  where
    N: Clone,
    E: Clone,
  {
    let (source_key, target_key) = {
      let edge = self
        .get_edge(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {} not found", edge_key))?;
      let edge = edge.read_arc();
      (edge.source(), edge.target())
    };

    let (target_inbound, target_outbound) = {
      let target_node = self
        .get_node(target_key)
        .ok_or_else(|| make_internal_report!("Target node {} not found", target_key))?;
      let target_node = target_node.read_arc();
      (target_node.inbound().to_vec(), target_node.outbound().to_vec())
    };

    for &inbound_edge_key in &target_inbound {
      if inbound_edge_key != edge_key {
        if let Some(inbound_edge) = self.get_edge(inbound_edge_key) {
          inbound_edge.write_arc().set_target(source_key);
          if let Some(source_node) = self.get_node(source_key) {
            let mut source_node = source_node.write_arc();
            if !source_node.inbound().contains(&inbound_edge_key) {
              source_node.inbound_mut().push(inbound_edge_key);
            }
          }
        }
      }
    }

    let mut new_edges = Vec::with_capacity(target_outbound.len());
    for &outbound_edge_key in &target_outbound {
      if outbound_edge_key != edge_key {
        if let Some(outbound_edge) = self.get_edge(outbound_edge_key) {
          new_edges.push(Arc::clone(&outbound_edge));
          outbound_edge.write_arc().set_source(source_key);
          if let Some(source_node) = self.get_node(source_key) {
            let mut source_node = source_node.write_arc();
            if !source_node.outbound().contains(&outbound_edge_key) {
              source_node.outbound_mut().push(outbound_edge_key);
            }
          }
        }
      }
    }

    if let Some(source_node) = self.get_node(source_key) {
      let mut source_node = source_node.write_arc();
      source_node.outbound_mut().retain(|&e| e != edge_key);
    }

    let removed_edge = self.remove_edge(edge_key)?;
    let (removed_node, _removed_edges) = self.remove_node(target_key)?;

    Ok((removed_node, removed_edge, new_edges))
  }
}

#[cfg(test)]
pub mod tests {
  use std::collections::BTreeMap;

  use pretty_assertions::assert_eq;

  use crate::graph::edge::Weighted;
  use crate::graph::node::Named;
  use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
  use crate::io::nwk::{
    EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight, nwk_read_str, nwk_write_str,
  };

  use super::*;

  pub type TestGraph = Graph<TestNode, TestEdge, ()>;

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
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
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
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
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
        .push(node.payload.name().unwrap().as_ref().to_owned());
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

  #[test]
  fn test_collapse_edge_simple_chain() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let internal_node = graph.add_node(TestNode(Some("internal".to_owned())));
    let leaf_node = graph.add_node(TestNode(Some("A".to_owned())));

    let root_to_internal_edge = graph.add_edge(root_node, internal_node, TestEdge(Some(0.1)))?;
    let internal_to_leaf_edge = graph.add_edge(internal_node, leaf_node, TestEdge(Some(0.2)))?;

    graph.build()?;

    graph.collapse_edge(root_to_internal_edge)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.2)root;");

    Ok(())
  }

  #[test]
  fn test_collapse_edge_binary_tree() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let internal_node = graph.add_node(TestNode(Some("internal".to_owned())));
    let a_node = graph.add_node(TestNode(Some("A".to_owned())));
    let b_node = graph.add_node(TestNode(Some("B".to_owned())));

    let root_to_internal_edge = graph.add_edge(root_node, internal_node, TestEdge(Some(0.1)))?;
    let internal_to_a_edge = graph.add_edge(internal_node, a_node, TestEdge(Some(0.2)))?;
    let root_to_b_edge = graph.add_edge(root_node, b_node, TestEdge(Some(0.3)))?;

    graph.build()?;

    graph.collapse_edge(root_to_internal_edge)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(B:0.3,A:0.2)root;");

    Ok(())
  }

  #[test]
  fn test_collapse_edge_complex_tree() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let left_internal = graph.add_node(TestNode(Some("left".to_owned())));
    let right_internal = graph.add_node(TestNode(Some("right".to_owned())));
    let a_node = graph.add_node(TestNode(Some("A".to_owned())));
    let b_node = graph.add_node(TestNode(Some("B".to_owned())));
    let c_node = graph.add_node(TestNode(Some("C".to_owned())));
    let d_node = graph.add_node(TestNode(Some("D".to_owned())));

    let root_to_left_edge = graph.add_edge(root_node, left_internal, TestEdge(Some(0.1)))?;
    let root_to_right_edge = graph.add_edge(root_node, right_internal, TestEdge(Some(0.2)))?;
    let left_to_a_edge = graph.add_edge(left_internal, a_node, TestEdge(Some(0.3)))?;
    let left_to_b_edge = graph.add_edge(left_internal, b_node, TestEdge(Some(0.4)))?;
    let right_to_c_edge = graph.add_edge(right_internal, c_node, TestEdge(Some(0.5)))?;
    let right_to_d_edge = graph.add_edge(right_internal, d_node, TestEdge(Some(0.6)))?;

    graph.build()?;

    graph.collapse_edge(root_to_left_edge)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "((C:0.5,D:0.6)right:0.2,A:0.3,B:0.4)root;");

    Ok(())
  }

  #[allow(clippy::assertions_on_result_states)]
  #[test]
  fn test_collapse_edge_invalid_edge() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let leaf_node = graph.add_node(TestNode(Some("A".to_owned())));
    graph.add_edge(root_node, leaf_node, TestEdge(Some(0.1)))?;
    graph.build()?;

    let invalid_key = GraphEdgeKey(9999);
    let result = graph.collapse_edge(invalid_key);

    assert!(result.is_err());

    Ok(())
  }

  #[test]
  fn test_collapse_edge_leaf_edge() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let a_node = graph.add_node(TestNode(Some("A".to_owned())));
    let b_node = graph.add_node(TestNode(Some("B".to_owned())));

    let root_to_a_edge = graph.add_edge(root_node, a_node, TestEdge(Some(0.1)))?;
    let root_to_b_edge = graph.add_edge(root_node, b_node, TestEdge(Some(0.2)))?;

    graph.build()?;

    graph.collapse_edge(root_to_a_edge)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(B:0.2)root;");

    Ok(())
  }

  #[test]
  fn test_collapse_edge_no_duplicate_edges() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let internal_node = graph.add_node(TestNode(Some("internal".to_owned())));
    let leaf_node = graph.add_node(TestNode(Some("A".to_owned())));

    // Create a scenario where source already has some connections that could duplicate
    let root_to_internal_edge = graph.add_edge(root_node, internal_node, TestEdge(Some(0.1)))?;
    let internal_to_leaf_edge = graph.add_edge(internal_node, leaf_node, TestEdge(Some(0.2)))?;

    graph.build()?;

    // Before collapse: root -> internal -> leaf
    let source_node_before = graph.get_node(root_node).unwrap();
    let initial_outbound_count = source_node_before.read_arc().outbound().len();

    graph.collapse_edge(root_to_internal_edge)?;

    // After collapse: root -> leaf (internal node removed)
    let source_node_after = graph.get_node(root_node).unwrap();
    let source_node_after = source_node_after.read_arc();

    // Verify no duplicate edges in adjacency lists
    let outbound_edges = source_node_after.outbound();
    let unique_outbound: HashSet<_> = outbound_edges.iter().collect();
    assert_eq!(
      outbound_edges.len(),
      unique_outbound.len(),
      "Duplicate outbound edges detected"
    );

    let inbound_edges = source_node_after.inbound();
    let unique_inbound: HashSet<_> = inbound_edges.iter().collect();
    assert_eq!(
      inbound_edges.len(),
      unique_inbound.len(),
      "Duplicate inbound edges detected"
    );

    // Should have one outbound edge (to leaf)
    assert_eq!(outbound_edges.len(), 1);

    Ok(())
  }

  #[test]
  fn test_collapse_edge_adjacency_lists_maintained() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let left_internal = graph.add_node(TestNode(Some("left".to_owned())));
    let right_internal = graph.add_node(TestNode(Some("right".to_owned())));
    let a_node = graph.add_node(TestNode(Some("A".to_owned())));
    let b_node = graph.add_node(TestNode(Some("B".to_owned())));

    let root_to_left_edge = graph.add_edge(root_node, left_internal, TestEdge(Some(0.1)))?;
    let root_to_right_edge = graph.add_edge(root_node, right_internal, TestEdge(Some(0.2)))?;
    let left_to_a_edge = graph.add_edge(left_internal, a_node, TestEdge(Some(0.3)))?;
    let left_to_b_edge = graph.add_edge(left_internal, b_node, TestEdge(Some(0.4)))?;

    graph.build()?;

    // Collapse root -> left_internal edge
    graph.collapse_edge(root_to_left_edge)?;

    // Verify root node's adjacency lists are correct
    let root_node_ref = graph.get_node(root_node).unwrap();
    let root_node_ref = root_node_ref.read_arc();

    // Root should have outbound edges to: right_internal, A, B
    assert_eq!(root_node_ref.outbound().len(), 3);

    // Verify all outbound edges from root point to correct targets
    let mut target_nodes = Vec::new();
    for &edge_key in root_node_ref.outbound() {
      if let Some(edge) = graph.get_edge(edge_key) {
        let target_key = edge.read_arc().target();
        if let Some(target_node) = graph.get_node(target_key) {
          if let Some(name) = target_node.read_arc().payload().read_arc().name() {
            target_nodes.push(name.as_ref().to_owned());
          }
        }
      }
    }
    target_nodes.sort();
    assert_eq!(target_nodes, vec!["A", "B", "right"]);

    // Verify leaf nodes have correct inbound edges
    let a_node_ref = graph.get_node(a_node).unwrap();
    let a_node_binding = a_node_ref.read_arc();
    let a_inbound = a_node_binding.inbound();
    assert_eq!(a_inbound.len(), 1);

    // Verify the edge from root to A has correct source
    if let Some(edge) = graph.get_edge(a_inbound[0]) {
      assert_eq!(edge.read_arc().source(), root_node);
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_multiple_inbound_edges() -> Result<(), Report> {
    // This test verifies handling when target node has multiple inbound edges
    // (though in a tree this shouldn't happen, we test the general case)
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let source1_node = graph.add_node(TestNode(Some("source1".to_owned())));
    let source2_node = graph.add_node(TestNode(Some("source2".to_owned())));
    let target_node = graph.add_node(TestNode(Some("target".to_owned())));
    let leaf_node = graph.add_node(TestNode(Some("leaf".to_owned())));

    // Create edges: source1 -> target, source2 -> target, target -> leaf
    let edge_to_collapse = graph.add_edge(source1_node, target_node, TestEdge(Some(0.1)))?;
    let other_inbound = graph.add_edge(source2_node, target_node, TestEdge(Some(0.2)))?;
    let outbound_edge = graph.add_edge(target_node, leaf_node, TestEdge(Some(0.3)))?;

    graph.build()?;

    // Collapse source1 -> target edge
    graph.collapse_edge(edge_to_collapse)?;

    // Verify source1 now has the outbound edge to leaf
    let source1_ref = graph.get_node(source1_node).unwrap();
    let source1_binding = source1_ref.read_arc();
    let source1_outbound = source1_binding.outbound();
    assert_eq!(source1_outbound.len(), 1);

    // Verify the edge now goes from source1 to leaf
    if let Some(edge) = graph.get_edge(source1_outbound[0]) {
      let edge_ref = edge.read_arc();
      assert_eq!(edge_ref.source(), source1_node);
      assert_eq!(edge_ref.target(), leaf_node);
    }

    // Verify source1 also inherited the other inbound edge (source2 -> source1)
    let source1_inbound = source1_binding.inbound();
    assert_eq!(source1_inbound.len(), 1);

    // Verify that source2 -> target edge now points to source1
    if let Some(edge) = graph.get_edge(source1_inbound[0]) {
      let edge_ref = edge.read_arc();
      assert_eq!(edge_ref.source(), source2_node);
      assert_eq!(edge_ref.target(), source1_node);
    }

    // Verify target node was removed
    assert!(graph.get_node(target_node).is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_edge_adjacency_consistency() -> Result<(), Report> {
    // Test that after edge collapse, all edge references in adjacency lists
    // correspond to actual edges that exist and point correctly
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let internal_node = graph.add_node(TestNode(Some("internal".to_owned())));
    let leaf1_node = graph.add_node(TestNode(Some("leaf1".to_owned())));
    let leaf2_node = graph.add_node(TestNode(Some("leaf2".to_owned())));

    let root_to_internal = graph.add_edge(root_node, internal_node, TestEdge(Some(0.1)))?;
    let internal_to_leaf1 = graph.add_edge(internal_node, leaf1_node, TestEdge(Some(0.2)))?;
    let internal_to_leaf2 = graph.add_edge(internal_node, leaf2_node, TestEdge(Some(0.3)))?;

    graph.build()?;

    graph.collapse_edge(root_to_internal)?;

    // Verify consistency: every edge key in node adjacency lists corresponds to a real edge
    for node in graph.nodes.iter().flatten() {
      let node_ref = node.read_arc();

      // Check outbound edges
      for &edge_key in node_ref.outbound() {
        let edge = graph.get_edge(edge_key).expect("Outbound edge should exist");
        let edge_ref = edge.read_arc();
        assert_eq!(edge_ref.source(), node_ref.key(), "Edge source should match node");
      }

      // Check inbound edges
      for &edge_key in node_ref.inbound() {
        let edge = graph.get_edge(edge_key).expect("Inbound edge should exist");
        let edge_ref = edge.read_arc();
        assert_eq!(edge_ref.target(), node_ref.key(), "Edge target should match node");
      }
    }

    Ok(())
  }
}
