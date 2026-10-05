use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use eyre::Report;
use itertools::Itertools;
use treetime_utils::{make_error, make_internal_report, make_report};

pub struct TreeView<'g> {
  graph: &'g Graph,
  root: GraphNodeKey,
  preorder: Vec<GraphNodeKey>,
  parents: Vec<Option<(GraphNodeKey, GraphEdgeKey)>>,
  child_offsets: Vec<usize>,
  children: Vec<(GraphNodeKey, GraphEdgeKey)>,
}

impl<'g> TreeView<'g> {
  pub fn new(graph: &'g Graph) -> Result<Self, Report> {
    let root = single_root(graph)?;
    let parents = single_parents(graph)?;
    let (child_offsets, children) = children_by_node(graph);
    let preorder = preorder_from(root, &child_offsets, &children);
    if preorder.len() != graph.get_nodes().count() {
      return Err(cycle_error(graph, &preorder, &parents));
    }
    Ok(Self {
      graph,
      root,
      preorder,
      parents,
      child_offsets,
      children,
    })
  }

  pub fn graph(&self) -> &'g Graph {
    self.graph
  }

  pub fn root(&self) -> GraphNodeKey {
    self.root
  }

  pub fn preorder(&self) -> &[GraphNodeKey] {
    &self.preorder
  }

  pub fn parent(&self, key: GraphNodeKey) -> Option<(GraphNodeKey, GraphEdgeKey)> {
    self.parents[key.as_usize()]
  }

  pub fn children(&self, key: GraphNodeKey) -> &[(GraphNodeKey, GraphEdgeKey)] {
    let index = key.as_usize();
    &self.children[self.child_offsets[index]..self.child_offsets[index + 1]]
  }
}

fn single_root(graph: &Graph) -> Result<GraphNodeKey, Report> {
  let roots = graph
    .get_nodes()
    .filter(|node| node.inbound().is_empty())
    .map(|node| node.key())
    .collect_vec();
  match roots.as_slice() {
    [root] => Ok(*root),
    [] => make_error!("The graph is not a tree: every node has a parent, so there is no root"),
    _ => make_error!(
      "The graph is not a tree: it has {} roots ({}), but a tree has one",
      roots.len(),
      roots.iter().map(|key| format!("node {key}")).join(", ")
    ),
  }
}

fn single_parents(graph: &Graph) -> Result<Vec<Option<(GraphNodeKey, GraphEdgeKey)>>, Report> {
  let mut parents = vec![None; graph.nodes.len()];
  for node in graph.get_nodes() {
    parents[node.key().as_usize()] = match node.inbound() {
      [] => None,
      [edge_key] => Some((graph.get_source_node_key(*edge_key)?, *edge_key)),
      inbound => {
        return make_error!(
          "The graph is not a tree: node {} has {} parents, but a node of a tree has at most one",
          node.key(),
          inbound.len()
        );
      },
    };
  }
  Ok(parents)
}

fn children_by_node(graph: &Graph) -> (Vec<usize>, Vec<(GraphNodeKey, GraphEdgeKey)>) {
  let mut offsets = Vec::with_capacity(graph.nodes.len() + 1);
  let mut children = Vec::with_capacity(graph.edges.len());
  offsets.push(0);
  for node in &graph.nodes {
    if let Some(node) = node {
      children.extend(graph.children_keys_of(node));
    }
    offsets.push(children.len());
  }
  (offsets, children)
}

fn preorder_from(
  root: GraphNodeKey,
  child_offsets: &[usize],
  children: &[(GraphNodeKey, GraphEdgeKey)],
) -> Vec<GraphNodeKey> {
  let mut order = Vec::with_capacity(child_offsets.len() - 1);
  let mut stack = vec![root];
  while let Some(key) = stack.pop() {
    order.push(key);
    let index = key.as_usize();
    let node_children = &children[child_offsets[index]..child_offsets[index + 1]];
    stack.extend(node_children.iter().rev().map(|(child_key, _)| *child_key));
  }
  order
}

fn cycle_error(graph: &Graph, preorder: &[GraphNodeKey], parents: &[Option<(GraphNodeKey, GraphEdgeKey)>]) -> Report {
  let mut reached = vec![false; parents.len()];
  for key in preorder {
    reached[key.as_usize()] = true;
  }
  let Some(mut key) = graph
    .get_nodes()
    .map(|node| node.key())
    .find(|key| !reached[key.as_usize()])
  else {
    return make_internal_report!("The tree preorder lists a node twice");
  };
  let mut seen = vec![false; parents.len()];
  while !seen[key.as_usize()] {
    seen[key.as_usize()] = true;
    let Some((parent_key, _)) = parents[key.as_usize()] else {
      return make_internal_report!("Node {key} has no parent, but it is not the root");
    };
    key = parent_key;
  }
  make_report!("The graph is not a tree: node {key} is on a cycle, which the root does not reach")
}
