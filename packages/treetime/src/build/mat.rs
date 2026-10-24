use crate::build::query::Query;
use crate::make_internal_report;
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Debug)]
pub struct Mat {
  graph: Graph,
  root: GraphNodeKey,
  reference: Seq,
  subs: Vec<Vec<Sub>>,
  samples: BTreeMap<GraphNodeKey, usize>,
  parsimony: usize,
  depth: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Detached {
  pub reattach_hint: [Option<GraphNodeKey>; 2],
}

impl Mat {
  pub fn new(reference: Seq) -> Self {
    let mut graph = Graph::new();
    let root = graph.add_node();
    Self {
      graph,
      root,
      reference,
      subs: vec![],
      samples: BTreeMap::new(),
      parsimony: 0,
      depth: vec![0],
    }
  }

  pub(crate) fn from_parts(
    graph: Graph,
    root: GraphNodeKey,
    reference: Seq,
    subs: Vec<Vec<Sub>>,
    samples: BTreeMap<GraphNodeKey, usize>,
  ) -> Self {
    let parsimony = graph
      .get_edges()
      .map(|edge| subs.get(edge.key().as_usize()).map_or(0, Vec::len))
      .sum();
    let mut mat = Self {
      graph,
      root,
      reference,
      subs,
      samples,
      parsimony,
      depth: vec![],
    };
    for node in mat.preorder().into_iter().rev() {
      let depth = mat
        .children(node)
        .map(|(child, edge)| mat.depth_below(child) + mat.edge_subs(edge).len())
        .max()
        .unwrap_or(0);
      mat.set_depth(node, depth);
    }
    mat
  }

  pub(crate) fn preorder(&self) -> Vec<GraphNodeKey> {
    let mut preorder = vec![];
    let mut stack = vec![self.root];
    while let Some(node) = stack.pop() {
      preorder.push(node);
      stack.extend(self.children(node).map(|(child, _)| child));
    }
    preorder
  }

  pub(crate) fn depth_below(&self, node: GraphNodeKey) -> usize {
    self.depth.get(node.as_usize()).copied().unwrap_or(0)
  }

  fn set_depth(&mut self, node: GraphNodeKey, depth: usize) {
    let index = node.as_usize();
    if self.depth.len() <= index {
      self.depth.resize(index + 1, 0);
    }
    self.depth[index] = depth;
  }

  fn raise_depth(&mut self, mut node: GraphNodeKey, mut depth: usize) {
    while depth > self.depth_below(node) {
      self.set_depth(node, depth);
      let Some((parent, edge)) = self.parent(node) else {
        break;
      };
      depth += self.edge_subs(edge).len();
      node = parent;
    }
  }

  pub fn graph(&self) -> &Graph {
    &self.graph
  }

  pub fn root(&self) -> GraphNodeKey {
    self.root
  }

  pub fn reference(&self) -> &Seq {
    &self.reference
  }

  pub fn parsimony(&self) -> usize {
    self.parsimony
  }

  pub fn samples(&self) -> &BTreeMap<GraphNodeKey, usize> {
    &self.samples
  }

  pub(crate) fn sample_of(&self, node: GraphNodeKey) -> Option<usize> {
    self.samples.get(&node).copied()
  }

  pub fn edge_subs(&self, edge: GraphEdgeKey) -> &[Sub] {
    self.subs.get(edge.as_usize()).map_or(&[], Vec::as_slice)
  }

  pub(crate) fn has_edge(&self, edge: GraphEdgeKey) -> bool {
    self.graph.get_edge(edge).is_some()
  }

  pub(crate) fn has_node(&self, node: GraphNodeKey) -> bool {
    self.graph.get_node(node).is_some()
  }

  pub(crate) fn is_leaf(&self, node: GraphNodeKey) -> bool {
    self.samples.contains_key(&node)
  }

  pub(crate) fn parent(&self, node: GraphNodeKey) -> Option<(GraphNodeKey, GraphEdgeKey)> {
    let &edge = self.graph.get_node(node)?.inbound().first()?;
    Some((self.graph[edge].source(), edge))
  }

  pub(crate) fn children(&self, node: GraphNodeKey) -> impl Iterator<Item = (GraphNodeKey, GraphEdgeKey)> + '_ {
    self
      .graph
      .get_node(node)
      .map(|node| node.outbound())
      .unwrap_or_default()
      .iter()
      .map(|&edge| (self.graph[edge].target(), edge))
  }

  pub(crate) fn n_children(&self, node: GraphNodeKey) -> usize {
    self.graph.get_node(node).map_or(0, |node| node.outbound().len())
  }

  pub(crate) fn edge_endpoints(&self, edge: GraphEdgeKey) -> (GraphNodeKey, GraphNodeKey) {
    let edge = &self.graph[edge];
    (edge.source(), edge.target())
  }

  pub(crate) fn node_state(&self, node: GraphNodeKey) -> Vec<(usize, AsciiChar)> {
    let mut path = vec![];
    let mut current = node;
    while let Some((parent, edge)) = self.parent(current) {
      path.push(edge);
      current = parent;
    }
    let mut state = vec![];
    for &edge in path.iter().rev() {
      self.apply_subs(&mut state, self.edge_subs(edge));
    }
    state
  }

  pub(crate) fn apply_subs<'a>(&self, state: &mut Vec<(usize, AsciiChar)>, subs: impl IntoIterator<Item = &'a Sub>) {
    for sub in subs {
      let index = state.partition_point(|&(pos, _)| pos < sub.pos());
      let present = state.get(index).is_some_and(|&(pos, _)| pos == sub.pos());
      let back_to_reference = sub.qry() == self.reference[sub.pos()];
      match (present, back_to_reference) {
        (true, true) => {
          state.remove(index);
        },
        (true, false) => state[index].1 = sub.qry(),
        (false, false) => state.insert(index, (sub.pos(), sub.qry())),
        (false, true) => {},
      }
    }
  }

  pub(crate) fn subs_to_query(&self, node_state: &[(usize, AsciiChar)], query: &Query) -> Result<Vec<Sub>, Report> {
    query
      .diffs()
      .iter()
      .merge_join_by(node_state, |a, b| a.0.cmp(&b.0))
      .filter_map(|entry| {
        let (pos, from, to) = match entry {
          itertools::EitherOrBoth::Both(&(pos, to), &(_, from)) => (pos, from, Some(to)),
          itertools::EitherOrBoth::Left(&(pos, to)) => (pos, self.reference[pos], Some(to)),
          itertools::EitherOrBoth::Right(&(pos, from)) => (pos, from, query.state(pos, &self.reference)),
        };
        let to = to?;
        (from != to).then(|| Sub::new(from, pos, to))
      })
      .collect()
  }

  fn set_subs(&mut self, edge: GraphEdgeKey, subs: Vec<Sub>) {
    let index = edge.as_usize();
    if self.subs.len() <= index {
      self.subs.resize_with(index + 1, Vec::new);
    }
    self.parsimony += subs.len();
    self.parsimony -= self.subs[index].len();
    self.subs[index] = subs;
  }

  fn take_subs(&mut self, edge: GraphEdgeKey) -> Vec<Sub> {
    let subs = self
      .subs
      .get_mut(edge.as_usize())
      .map(std::mem::take)
      .unwrap_or_default();
    self.parsimony -= subs.len();
    subs
  }

  fn connect(&mut self, parent: GraphNodeKey, child: GraphNodeKey, subs: Vec<Sub>) -> Result<GraphEdgeKey, Report> {
    let edge = self.graph.add_edge(parent, child)?;
    let depth = self.depth_below(child) + subs.len();
    self.set_subs(edge, subs);
    self.raise_depth(parent, depth);
    Ok(edge)
  }

  pub(crate) fn add_sample(
    &mut self,
    parent: GraphNodeKey,
    subs: Vec<Sub>,
    sample: usize,
  ) -> Result<GraphNodeKey, Report> {
    let leaf = self.graph.add_node();
    self.connect(parent, leaf, subs)?;
    self.samples.insert(leaf, sample);
    Ok(leaf)
  }

  pub(crate) fn attach(&mut self, parent: GraphNodeKey, child: GraphNodeKey, subs: Vec<Sub>) -> Result<(), Report> {
    self.connect(parent, child, subs)?;
    Ok(())
  }

  pub(crate) fn split_edge(&mut self, edge: GraphEdgeKey, upper: &[usize]) -> Result<GraphNodeKey, Report> {
    let (parent, child) = self.edge_endpoints(edge);
    let (upper_subs, lower_subs): (Vec<Sub>, Vec<Sub>) = self
      .take_subs(edge)
      .into_iter()
      .partition(|sub| upper.binary_search(&sub.pos()).is_ok());
    self.graph.remove_edge(edge)?;
    let node = self.graph.add_node();
    self.connect(parent, node, upper_subs)?;
    self.connect(node, child, lower_subs)?;
    Ok(node)
  }

  pub(crate) fn detach(&mut self, subtree: GraphNodeKey) -> Result<Detached, Report> {
    let (parent, edge) = self
      .parent(subtree)
      .ok_or_else(|| make_internal_report!("Cannot detach the root {subtree}"))?;
    self.take_subs(edge);
    self.graph.remove_edge(edge)?;
    let reattach_hint = match self.collapse_unary(parent)? {
      Some((grandparent, child)) => [Some(grandparent), Some(child)],
      None => [Some(parent), None],
    };
    Ok(Detached { reattach_hint })
  }

  fn collapse_unary(&mut self, node: GraphNodeKey) -> Result<Option<(GraphNodeKey, GraphNodeKey)>, Report> {
    if self.is_leaf(node) || self.n_children(node) != 1 {
      return Ok(None);
    }
    let Some((grandparent, upper_edge)) = self.parent(node) else {
      return Ok(None);
    };
    let Some((child, lower_edge)) = self.children(node).next() else {
      return Ok(None);
    };
    let upper = self.take_subs(upper_edge);
    let lower = self.take_subs(lower_edge);
    let composed = compose_substitutions(&upper, &lower)?;
    self.graph.remove_node(node)?;
    self.connect(grandparent, child, composed)?;
    Ok(Some((grandparent, child)))
  }

  pub(crate) fn collapse_empty_internal_edges(&mut self) -> Result<usize, Report> {
    let empty = self
      .graph
      .get_edges()
      .filter(|edge| self.edge_subs(edge.key()).is_empty() && !self.is_leaf(edge.target()))
      .map(|edge| edge.key())
      .collect_vec();
    let mut collapsed = 0;
    for edge in empty {
      if !self.has_edge(edge) {
        continue;
      }
      let (parent, child) = self.edge_endpoints(edge);
      let grandchildren = self.children(child).map(|(_, edge)| edge).collect_vec();
      for grandchild_edge in grandchildren {
        self.graph.reparent_edge(grandchild_edge, parent)?;
      }
      self.graph.remove_node(child)?;
      collapsed += 1;
    }
    Ok(collapsed)
  }

  pub(crate) fn into_parts(self) -> (Graph, Vec<Vec<Sub>>, BTreeMap<GraphNodeKey, usize>) {
    (self.graph, self.subs, self.samples)
  }
}
