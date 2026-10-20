use crate::edge::GraphEdgeKey;
use crate::graph::{Graph, missing_key};
use crate::node::{GraphNodeKey, Node};
use eyre::Report;
use itertools::Itertools;
use std::mem;
use treetime_utils::{make_error, make_internal_report, make_report};

pub struct GraphSchedule<'g> {
  graph: &'g Graph,
  roots: Vec<GraphNodeKey>,
  leaves: Vec<GraphNodeKey>,
  preorder: Vec<GraphNodeKey>,
  postorder: Vec<GraphNodeKey>,
  parents: Adjacency,
  children: Adjacency,
}

impl<'g> GraphSchedule<'g> {
  pub fn new(graph: &'g Graph) -> Result<Self, Report> {
    let parents = Adjacency::build(graph, |node| graph.parents_keys_of(node));
    let children = Adjacency::build(graph, |node| graph.children_keys_of(node));
    let roots = roots_of(graph, &parents);
    if roots.is_empty() {
      return make_error!("The graph has no root: every node has a parent");
    }
    let preorder = preorder_of(graph, &roots, &parents, &children);
    if preorder.len() != graph.get_nodes().count() {
      return Err(cycle_error(graph, &preorder, &parents));
    }
    let postorder = postorder_of(graph, &roots, &children);
    let leaves = preorder
      .iter()
      .copied()
      .filter(|&key| children.of(key).is_empty())
      .collect_vec();
    Ok(Self {
      graph,
      roots,
      leaves,
      preorder,
      postorder,
      parents,
      children,
    })
  }

  pub fn preorder(&self) -> &[GraphNodeKey] {
    &self.preorder
  }

  pub fn postorder(&self) -> &[GraphNodeKey] {
    &self.postorder
  }

  pub fn parents(&self, key: GraphNodeKey) -> &[(GraphNodeKey, GraphEdgeKey)] {
    self.parents.of(key)
  }

  pub fn children(&self, key: GraphNodeKey) -> &[(GraphNodeKey, GraphEdgeKey)] {
    self.children.of(key)
  }

  pub fn roots(&self) -> &[GraphNodeKey] {
    &self.roots
  }

  pub fn leaves(&self) -> &[GraphNodeKey] {
    &self.leaves
  }

  pub fn child_slot(&self, edge: GraphEdgeKey) -> usize {
    self.children.slot(edge)
  }

  pub fn parent_slot(&self, edge: GraphEdgeKey) -> usize {
    self.parents.slot(edge)
  }

  pub fn require_tree(&self) -> Result<(), Report> {
    if self.roots.len() != 1 {
      return make_error!(
        "The graph is not a tree: it has {} roots ({}), but a tree has one",
        self.roots.len(),
        self.roots.iter().map(|key| format!("node {key}")).join(", ")
      );
    }
    if let Some(key) = self.graph.node_keys().find(|&key| self.parents(key).len() > 1) {
      return make_error!(
        "The graph is not a tree: node {key} has {} parents, but a node of a tree has at most one",
        self.parents(key).len()
      );
    }
    Ok(())
  }

  pub fn into_tree(self) -> Result<TreeSchedule<'g>, Report> {
    self.require_tree()?;
    Ok(TreeSchedule { schedule: self })
  }

  pub fn contains(&self, key: GraphNodeKey) -> bool {
    self.graph.get_node(key).is_some()
  }

  pub(crate) fn contains_edge(&self, edge: GraphEdgeKey) -> bool {
    self.parents.slots.get(edge.as_usize()).is_some_and(Option::is_some)
  }

  pub(crate) fn edge_count(&self) -> usize {
    self.parents.entries.len()
  }

  pub(crate) fn key_capacity(&self) -> usize {
    self.graph.nodes.len()
  }
}

pub struct TreeSchedule<'g> {
  schedule: GraphSchedule<'g>,
}

impl<'g> TreeSchedule<'g> {
  pub fn root(&self) -> GraphNodeKey {
    self.schedule.roots[0]
  }

  pub fn parent(&self, key: GraphNodeKey) -> Option<(GraphNodeKey, GraphEdgeKey)> {
    self.schedule.parents(key).first().copied()
  }

  pub fn children(&self, key: GraphNodeKey) -> &[(GraphNodeKey, GraphEdgeKey)] {
    self.schedule.children(key)
  }

  pub fn preorder(&self) -> &[GraphNodeKey] {
    self.schedule.preorder()
  }

  pub fn schedule(&self) -> &GraphSchedule<'g> {
    &self.schedule
  }
}

struct Adjacency {
  offsets: Vec<usize>,
  entries: Vec<(GraphNodeKey, GraphEdgeKey)>,
  slots: Vec<Option<usize>>,
}

impl Adjacency {
  fn build<'a, I>(graph: &'a Graph, neighbors: impl Fn(&'a Node) -> I) -> Self
  where
    I: Iterator<Item = (GraphNodeKey, GraphEdgeKey)>,
  {
    let mut offsets = Vec::with_capacity(graph.nodes.len() + 1);
    let mut entries = Vec::with_capacity(graph.edges.len());
    let mut slots = vec![None; graph.edges.len()];
    offsets.push(0);
    for node in &graph.nodes {
      if let Some(node) = node {
        for (slot, (key, edge_key)) in neighbors(node).enumerate() {
          slots[edge_key.as_usize()] = Some(slot);
          entries.push((key, edge_key));
        }
      }
      offsets.push(entries.len());
    }
    Self {
      offsets,
      entries,
      slots,
    }
  }

  fn of(&self, key: GraphNodeKey) -> &[(GraphNodeKey, GraphEdgeKey)] {
    let index = key.as_usize();
    &self.entries[self.offsets[index]..self.offsets[index + 1]]
  }

  fn slot(&self, edge: GraphEdgeKey) -> usize {
    self
      .slots
      .get(edge.as_usize())
      .copied()
      .flatten()
      .unwrap_or_else(|| missing_key("Edge", edge))
  }
}

fn roots_of(graph: &Graph, parents: &Adjacency) -> Vec<GraphNodeKey> {
  let is_root = |key: GraphNodeKey| graph.get_node(key).is_some() && parents.of(key).is_empty();
  let mut listed = vec![false; graph.nodes.len()];
  let mut roots = Vec::with_capacity(graph.roots.len());
  for &key in &graph.roots {
    if is_root(key) && !mem::replace(&mut listed[key.as_usize()], true) {
      roots.push(key);
    }
  }
  roots.extend(graph.node_keys().filter(|&key| is_root(key) && !listed[key.as_usize()]));
  roots
}

fn preorder_of(graph: &Graph, roots: &[GraphNodeKey], parents: &Adjacency, children: &Adjacency) -> Vec<GraphNodeKey> {
  let mut waiting = (0..graph.nodes.len())
    .map(|index| parents.of(GraphNodeKey(index)).len())
    .collect_vec();
  let mut order = Vec::with_capacity(graph.nodes.len());
  let mut stack = roots.iter().rev().copied().collect_vec();
  while let Some(key) = stack.pop() {
    order.push(key);
    for &(child, _) in children.of(key).iter().rev() {
      let count = &mut waiting[child.as_usize()];
      *count -= 1;
      if *count == 0 {
        stack.push(child);
      }
    }
  }
  order
}

fn cycle_error(graph: &Graph, preorder: &[GraphNodeKey], parents: &Adjacency) -> Report {
  let mut visited = vec![false; graph.nodes.len()];
  for key in preorder {
    visited[key.as_usize()] = true;
  }
  let Some(mut key) = graph.node_keys().find(|key| !visited[key.as_usize()]) else {
    return make_internal_report!("The graph preorder lists a node twice");
  };
  let mut seen = vec![false; graph.nodes.len()];
  while !mem::replace(&mut seen[key.as_usize()], true) {
    let Some(&(parent, _)) = parents.of(key).iter().find(|(parent, _)| !visited[parent.as_usize()]) else {
      return make_internal_report!("Node {key} is not in the preorder, but every parent of it is");
    };
    key = parent;
  }
  make_report!("The graph has a cycle through node {key}")
}

fn postorder_of(graph: &Graph, roots: &[GraphNodeKey], children: &Adjacency) -> Vec<GraphNodeKey> {
  let mut entered = vec![false; graph.nodes.len()];
  let mut order = Vec::with_capacity(graph.nodes.len());
  let mut stack = roots.iter().rev().map(|&key| Frame::Enter(key)).collect_vec();
  while let Some(frame) = stack.pop() {
    match frame {
      Frame::Enter(key) => {
        if mem::replace(&mut entered[key.as_usize()], true) {
          continue;
        }
        stack.push(Frame::Exit(key));
        stack.extend(children.of(key).iter().rev().map(|&(child, _)| Frame::Enter(child)));
      },
      Frame::Exit(key) => order.push(key),
    }
  }
  order
}

enum Frame {
  Enter(GraphNodeKey),
  Exit(GraphNodeKey),
}
