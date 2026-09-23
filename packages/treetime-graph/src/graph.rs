use crate::edge::{Edge, GraphEdgeKey};
use crate::node::{GraphNodeKey, Node};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use treetime_utils::{make_internal_error, make_internal_report};

#[expect(
  clippy::field_scoped_visibility_modifiers,
  reason = "graph storage stays crate-internal behind accessor methods"
)]
#[derive(Debug, Serialize, Deserialize)]
pub struct Graph {
  pub(crate) nodes: Vec<Option<Node>>,
  pub(crate) edges: Vec<Option<Edge>>,
  pub(crate) roots: Vec<GraphNodeKey>,
  pub(crate) leaves: Vec<GraphNodeKey>,
}

impl Graph {
  pub fn new() -> Self {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      roots: vec![],
      leaves: vec![],
    }
  }

  pub fn get_node(&self, key: GraphNodeKey) -> Option<&Node> {
    self.nodes.get(key.as_usize())?.as_ref()
  }

  pub(crate) fn get_node_mut(&mut self, key: GraphNodeKey) -> Option<&mut Node> {
    self.nodes.get_mut(key.as_usize())?.as_mut()
  }

  pub fn get_edge(&self, key: GraphEdgeKey) -> Option<&Edge> {
    self.edges.get(key.as_usize())?.as_ref()
  }

  pub(crate) fn get_edge_mut(&mut self, key: GraphEdgeKey) -> Option<&mut Edge> {
    self.edges.get_mut(key.as_usize())?.as_mut()
  }

  pub fn parents_of<'a>(&'a self, node: &'a Node) -> impl DoubleEndedIterator<Item = (&'a Node, &'a Edge)> + 'a {
    node.inbound().iter().filter_map(move |&edge_key| {
      let edge = self.get_edge(edge_key)?;
      Some((self.get_node(edge.source())?, edge))
    })
  }

  pub(crate) fn parents_keys_of<'a>(
    &'a self,
    node: &'a Node,
  ) -> impl DoubleEndedIterator<Item = (GraphNodeKey, GraphEdgeKey)> + 'a {
    node
      .inbound()
      .iter()
      .filter_map(move |&edge_key| self.get_edge(edge_key).map(|edge| (edge.source(), edge_key)))
  }

  pub fn exactly_one_parent_of(&self, node: &Node) -> Result<(GraphNodeKey, GraphEdgeKey), Report> {
    self.one_parent_of(node)?.ok_or_else(|| {
      make_internal_report!(
        "No parents found for node {} (context: is_root={} is_leaf={})",
        node.key(),
        node.is_root(),
        node.is_leaf()
      )
    })
  }

  fn one_parent_of(&self, node: &Node) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report> {
    match node.inbound().len() {
      0 => Ok(None),
      1 => {
        let edge_key = node.inbound()[0];
        Ok(Some((self.get_source_node_key(edge_key)?, edge_key)))
      },
      n => make_internal_error!(
        "Only trees with exactly one parent per node are currently supported, but node '{}' has {n} parents",
        node.key()
      ),
    }
  }

  pub fn children_of<'a>(&'a self, node: &'a Node) -> impl DoubleEndedIterator<Item = (&'a Node, &'a Edge)> + 'a {
    node.outbound().iter().filter_map(move |&edge_key| {
      let edge = self.get_edge(edge_key)?;
      Some((self.get_node(edge.target())?, edge))
    })
  }

  pub fn children_keys_of<'a>(
    &'a self,
    node: &'a Node,
  ) -> impl DoubleEndedIterator<Item = (GraphNodeKey, GraphEdgeKey)> + 'a {
    node
      .outbound()
      .iter()
      .filter_map(move |&edge_key| self.get_edge(edge_key).map(|edge| (edge.target(), edge_key)))
  }

  pub fn get_source_node_key(&self, edge_key: GraphEdgeKey) -> Result<GraphNodeKey, Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    Ok(edge.source())
  }

  pub fn get_target_node_key(&self, edge_key: GraphEdgeKey) -> Result<GraphNodeKey, Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    Ok(edge.target())
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

  pub fn get_nodes(&self) -> impl DoubleEndedIterator<Item = &Node> {
    self.nodes.iter().filter_map(Option::as_ref)
  }

  pub(crate) fn node_keys(&self) -> impl DoubleEndedIterator<Item = GraphNodeKey> + '_ {
    self.get_nodes().map(Node::key)
  }

  pub fn get_exactly_one_root(&self) -> Result<&Node, Report> {
    if self.roots.len() != 1 {
      return make_internal_error!(
        "Only trees with exactly one root are currently supported, but found {} roots",
        self.roots.len()
      );
    }
    self
      .get_node(self.roots[0])
      .ok_or_else(|| make_internal_report!("Root node {} not found", self.roots[0]))
  }

  pub fn get_roots(&self) -> impl Iterator<Item = &Node> + '_ {
    self.roots.iter().filter_map(|key| self.get_node(*key))
  }

  pub fn root_keys(&self) -> impl Iterator<Item = GraphNodeKey> + '_ {
    self.roots.iter().copied()
  }

  pub fn get_leaves(&self) -> impl Iterator<Item = &Node> + '_ {
    self.leaves.iter().filter_map(|key| self.get_node(*key))
  }

  pub fn leaf_keys(&self) -> impl Iterator<Item = GraphNodeKey> + '_ {
    self.leaves.iter().copied()
  }

  pub fn get_internal_nodes(&self) -> impl DoubleEndedIterator<Item = &Node> {
    self.get_nodes().filter(|node| !node.is_leaf())
  }

  pub fn get_inner_nodes(&self) -> impl DoubleEndedIterator<Item = &Node> {
    self.get_nodes().filter(|node| !node.is_leaf() && !node.is_root())
  }

  pub fn get_edges(&self) -> impl DoubleEndedIterator<Item = &Edge> {
    self.edges.iter().filter_map(Option::as_ref)
  }

  pub fn edge_keys(&self) -> impl DoubleEndedIterator<Item = GraphEdgeKey> + '_ {
    self.get_edges().map(Edge::key)
  }

  pub(crate) fn path_from_root_to_node(&self, node_key: GraphNodeKey) -> Result<Vec<GraphNodeKey>, Report> {
    let mut node = self
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node not found on the graph: {node_key}"))?;

    let mut path = vec![node.key()];
    while let Some((parent_key, _)) = self.one_parent_of(node)? {
      path.push(parent_key);
      node = self
        .get_node(parent_key)
        .ok_or_else(|| make_internal_report!("Parent node not found on the graph: {parent_key}"))?;
    }

    path.reverse();
    Ok(path)
  }

  pub(crate) fn path_from_node_to_node(
    &self,
    start: GraphNodeKey,
    finish: GraphNodeKey,
  ) -> Result<Vec<(GraphNodeKey, Option<GraphEdgeKey>)>, Report> {
    let mut node = self
      .get_node(start)
      .ok_or_else(|| make_internal_report!("Node not found on the graph: {start}"))?;

    let mut path = vec![(node.key(), None)];
    loop {
      match self.one_parent_of(node)? {
        None => {
          return make_internal_error!(
            "When searching path from starting node {start} to destination node {finish}: reached root node without finding the destination"
          );
        },

        Some((parent_key, edge_key)) => {
          path.push((node.key(), Some(edge_key)));

          if parent_key == finish {
            break;
          }

          node = self
            .get_node(parent_key)
            .ok_or_else(|| make_internal_report!("Parent node not found on the graph: {parent_key}"))?;
        },
      }
    }

    Ok(path)
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

  pub fn degree_out(&self, key: GraphNodeKey) -> Result<usize, Report> {
    self
      .get_node(key)
      .map(Node::degree_out)
      .ok_or_else(|| make_internal_report!("Node not found: {key}"))
  }

  pub(crate) fn degree_in(&self, key: GraphNodeKey) -> Result<usize, Report> {
    self
      .get_node(key)
      .map(Node::degree_in)
      .ok_or_else(|| make_internal_report!("Node not found: {key}"))
  }

  pub fn has_parents(&self, node_key: GraphNodeKey) -> bool {
    self.get_node(node_key).is_some_and(Node::has_parents)
  }

  pub fn has_one_parent(&self, node_key: GraphNodeKey) -> bool {
    self.get_node(node_key).is_some_and(Node::has_one_parent)
  }

  pub fn has_at_most_one_parent(&self, node_key: GraphNodeKey) -> bool {
    self.get_node(node_key).is_some_and(Node::has_at_most_one_parent)
  }

  pub fn has_children(&self, node_key: GraphNodeKey) -> bool {
    self.get_node(node_key).is_some_and(Node::has_children)
  }

  pub fn has_one_child(&self, node_key: GraphNodeKey) -> bool {
    self.get_node(node_key).is_some_and(Node::has_one_child)
  }

  pub fn has_at_most_one_child(&self, node_key: GraphNodeKey) -> bool {
    self.get_node(node_key).is_some_and(Node::has_at_most_one_child)
  }

  pub fn parent_inbound_edge(&self, key: GraphNodeKey) -> Result<Option<GraphEdgeKey>, Report> {
    let node = self
      .get_node(key)
      .ok_or_else(|| make_internal_report!("Node not found: {key}"))?;
    Ok(node.inbound().first().copied())
  }

  pub fn edge_endpoints(&self, edge_key: GraphEdgeKey) -> Result<(GraphNodeKey, GraphNodeKey), Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    Ok((edge.source(), edge.target()))
  }

  pub fn node_parent(&self, node_key: GraphNodeKey) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report> {
    let node = self
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found"))?;
    self.one_parent_of(node)
  }

  pub fn root_key(&self) -> Result<GraphNodeKey, Report> {
    Ok(self.get_exactly_one_root()?.key())
  }
}

impl Default for Graph {
  fn default() -> Self {
    Self::new()
  }
}
