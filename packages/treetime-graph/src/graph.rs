use crate::edge::{Edge, GraphEdge, GraphEdgeKey};
use crate::graph_traverse::GraphNodeSafe;
use crate::node::{GraphNode, GraphNodeKey, Node};
use eyre::Report;
use itertools::Itertools;
use parking_lot::lock_api::{ArcRwLockReadGuard, ArcRwLockWriteGuard};
use parking_lot::{RawRwLock, RwLock};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::sync::Arc;
use treetime_utils::{make_internal_error, make_internal_report};

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

#[allow(clippy::field_scoped_visibility_modifiers)]
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph<N, E, D = ()>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  pub(crate) nodes: Vec<Option<Arc<RwLock<Node<N>>>>>,
  pub(crate) edges: Vec<Option<Arc<RwLock<Edge<E>>>>>,
  pub(crate) roots: Vec<GraphNodeKey>,
  pub(crate) leaves: Vec<GraphNodeKey>,
  pub(crate) data: Arc<RwLock<D>>,
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
        parents.len()
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
}
