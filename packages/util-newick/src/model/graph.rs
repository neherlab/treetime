use crate::model::data::{NewickEdgeData, NewickNodeData};
use crate::model::equality::graphs_equal;
use crate::model::traverse::{Postorder, Preorder};
use crate::model::validate::{describe_node, validate_graph};
use eyre::{Report, eyre};
use deser::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickGraph {
  nodes: Vec<NodeEntry>,
  edges: Vec<NewickEdgeEntry>,
  root: usize,
  rooted: Option<bool>,
  weight: Option<f64>,
  root_edge: NewickEdgeData,
}

impl NewickGraph {
  pub fn new(root: NewickNodeData) -> Self {
    Self {
      nodes: vec![NodeEntry::new(root)],
      edges: Vec::new(),
      root: 0,
      rooted: None,
      weight: None,
      root_edge: NewickEdgeData::new(),
    }
  }

  pub fn add_node(&mut self, data: NewickNodeData) -> usize {
    self.nodes.push(NodeEntry::new(data));
    self.nodes.len() - 1
  }

  pub fn add_edge(&mut self, parent: usize, child: usize, data: NewickEdgeData) -> Result<usize, Report> {
    let node_count = self.nodes.len();
    if parent >= node_count || child >= node_count {
      return Err(eyre!(
        "Cannot connect node {parent} to node {child}: the graph has {node_count} nodes"
      ));
    }
    if child == self.root {
      return Err(eyre!("The root, {}, cannot have a parent", describe_node(self, child)));
    }
    if parent == child {
      return Err(eyre!("{} cannot be its own parent", describe_node(self, child)));
    }
    if self.children(parent).any(|existing| existing == child) {
      return Err(eyre!(
        "{} already has an edge to {}",
        describe_node(self, parent),
        describe_node(self, child)
      ));
    }
    let child_entry = &self.nodes[child];
    if child_entry.data.hybrid().is_none() && !child_entry.parents.is_empty() {
      return Err(eyre!(
        "{} already has a parent, and only a hybrid node can have more than one",
        describe_node(self, child)
      ));
    }
    Ok(self.push_edge(parent, child, data))
  }

  pub fn add_child(&mut self, parent: usize, edge: NewickEdgeData, node: NewickNodeData) -> Result<usize, Report> {
    if parent >= self.nodes.len() {
      return Err(eyre!(
        "Cannot add a child to node {parent}: the graph has {} nodes",
        self.nodes.len()
      ));
    }
    let child = self.add_node(node);
    self.push_edge(parent, child, edge);
    Ok(child)
  }

  pub fn root(&self) -> usize {
    self.root
  }

  pub fn node_count(&self) -> usize {
    self.nodes.len()
  }

  pub fn edge_count(&self) -> usize {
    self.edges.len()
  }

  pub fn node(&self, node: usize) -> &NewickNodeData {
    &self.nodes[node].data
  }

  pub fn node_mut(&mut self, node: usize) -> &mut NewickNodeData {
    &mut self.nodes[node].data
  }

  pub fn nodes(&self) -> impl Iterator<Item = (usize, &NewickNodeData)> {
    self.nodes.iter().enumerate().map(|(idx, entry)| (idx, &entry.data))
  }

  pub fn edge(&self, edge: usize) -> &NewickEdgeEntry {
    &self.edges[edge]
  }

  pub fn edge_mut(&mut self, edge: usize) -> &mut NewickEdgeData {
    &mut self.edges[edge].data
  }

  pub fn edges(&self) -> impl Iterator<Item = (usize, &NewickEdgeEntry)> {
    self.edges.iter().enumerate()
  }

  pub fn child_edges(&self, node: usize) -> &[usize] {
    &self.nodes[node].children
  }

  pub fn parent_edges(&self, node: usize) -> &[usize] {
    &self.nodes[node].parents
  }

  pub fn children(&self, node: usize) -> impl Iterator<Item = usize> {
    self.nodes[node].children.iter().map(|&edge| self.edges[edge].child)
  }

  pub fn parents(&self, node: usize) -> impl Iterator<Item = usize> {
    self.nodes[node].parents.iter().map(|&edge| self.edges[edge].parent)
  }

  pub fn is_leaf(&self, node: usize) -> bool {
    self.nodes[node].children.is_empty()
  }

  pub fn rooted(&self) -> Option<bool> {
    self.rooted
  }

  pub fn set_rooted(&mut self, rooted: Option<bool>) {
    self.rooted = rooted;
  }

  pub fn weight(&self) -> Option<f64> {
    self.weight
  }

  pub fn set_weight(&mut self, weight: Option<f64>) {
    self.weight = weight;
  }

  pub fn root_edge(&self) -> &NewickEdgeData {
    &self.root_edge
  }

  pub fn root_edge_mut(&mut self) -> &mut NewickEdgeData {
    &mut self.root_edge
  }

  pub fn preorder(&self) -> Preorder<'_> {
    Preorder::new(self)
  }

  pub fn postorder(&self) -> Postorder<'_> {
    Postorder::new(self)
  }

  pub fn validate(&self) -> Result<(), Report> {
    validate_graph(self)
  }

  pub fn eq_ordered(&self, other: &Self) -> bool {
    graphs_equal(self, other, true)
  }

  pub(crate) fn from_parts(nodes: Vec<NewickNodeData>, edges: Vec<NewickEdgeEntry>, root: usize) -> Self {
    let mut entries: Vec<NodeEntry> = nodes.into_iter().map(NodeEntry::new).collect();
    for (idx, edge) in edges.iter().enumerate() {
      entries[edge.parent].children.push(idx);
      entries[edge.child].parents.push(idx);
    }
    Self {
      nodes: entries,
      edges,
      root,
      rooted: None,
      weight: None,
      root_edge: NewickEdgeData::new(),
    }
  }

  pub(crate) fn set_root_edge(&mut self, edge: NewickEdgeData) {
    self.root_edge = edge;
  }

  fn push_edge(&mut self, parent: usize, child: usize, data: NewickEdgeData) -> usize {
    let idx = self.edges.len();
    self.edges.push(NewickEdgeEntry { parent, child, data });
    self.nodes[parent].children.push(idx);
    self.nodes[child].parents.push(idx);
    idx
  }
}

impl PartialEq for NewickGraph {
  fn eq(&self, other: &Self) -> bool {
    graphs_equal(self, other, false)
  }
}

impl Eq for NewickGraph {}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickEdgeEntry {
  parent: usize,
  child: usize,
  data: NewickEdgeData,
}

impl NewickEdgeEntry {
  pub(crate) fn new(parent: usize, child: usize, data: NewickEdgeData) -> Self {
    Self { parent, child, data }
  }

  pub fn parent(&self) -> usize {
    self.parent
  }

  pub fn child(&self) -> usize {
    self.child
  }

  pub fn data(&self) -> &NewickEdgeData {
    &self.data
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct NodeEntry {
  data: NewickNodeData,
  children: Vec<usize>,
  parents: Vec<usize>,
}

impl NodeEntry {
  fn new(data: NewickNodeData) -> Self {
    Self {
      data,
      children: Vec::new(),
      parents: Vec::new(),
    }
  }
}
