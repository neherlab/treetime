use crate::breadth_first::{
  GraphTraversalContinuation, directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward,
};
use crate::edge::{GraphEdge, GraphEdgeKey};
use crate::graph::{Graph, NodeEdgePair, NodeEdgePayloadPair, SafeEdgePayloadRefMut, SafeNodePayloadRefMut};
use crate::node::{GraphNode, GraphNodeKey, Node};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::RwLock;
use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;
use traversal::Bft;
use treetime_utils::collections::container::{get_exactly_one, get_exactly_one_mut};

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

impl<N, E, D> Graph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
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

  /// Synchronously traverse graph in depth-first postorder fashion forward (from roots to leaves, along edge directions).
  ///
  /// Guarantees that for each visited node, all of it children (recursively) are visited before
  /// the node itself is visited.
  pub fn iter_depth_first_postorder_forward(&self, mut explorer: impl FnMut(GraphNodeBackward<N, E, D>)) {
    let root = self.get_exactly_one_root().unwrap();
    let mut stack = Vec::new();
    let mut visited = BTreeSet::new();
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
}
