use crate::equality::graphs_equal;
use crate::validate::validate_graph;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, Default, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickReadOptions {
  pub enewick: bool,
}

#[derive(Clone, Debug, SmartDefault, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickWriteOptions {
  pub style: NwkStyle,
  pub significant_digits: Option<u8>,
  pub decimal_digits: Option<i8>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, SmartDefault, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum NwkStyle {
  Plain,
  #[default]
  Beast,
  Nhx,
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NexusTree {
  pub name: String,
  pub graph: NewickGraph,
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickGraph {
  pub nodes: Vec<NewickNodeData>,
  pub edges: Vec<NewickEdgeEntry>,
  pub root: usize,
  pub rooted: Option<bool>,
}

impl NewickGraph {
  pub fn new() -> Self {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      root: 0,
      rooted: None,
    }
  }

  pub fn add_node(&mut self, data: NewickNodeData) -> usize {
    let idx = self.nodes.len();
    self.nodes.push(data);
    idx
  }

  pub fn add_edge(&mut self, parent: usize, child: usize, data: NewickEdgeData) -> usize {
    let idx = self.edges.len();
    self.edges.push(NewickEdgeEntry { parent, child, data });
    self.nodes[parent].children.push(idx);
    idx
  }

  pub fn validate(&self) -> Result<(), Report> {
    validate_graph(self)
  }

  pub fn eq_ordered(&self, other: &NewickGraph) -> bool {
    graphs_equal(self, other, true)
  }
}

impl Default for NewickGraph {
  fn default() -> Self {
    Self::new()
  }
}

impl PartialEq for NewickGraph {
  fn eq(&self, other: &Self) -> bool {
    graphs_equal(self, other, false)
  }
}

impl Eq for NewickGraph {}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickEdgeEntry {
  pub parent: usize,
  pub child: usize,
  pub data: NewickEdgeData,
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickNodeData {
  pub label: Option<NewickLabel>,
  pub node_attrs: BTreeMap<String, NewickValue>,
  pub raw_comments: Vec<String>,
  pub hybrid: Option<NewickHybrid>,
  pub children: Vec<usize>,
}

impl NewickNodeData {
  pub fn new() -> Self {
    Self {
      label: None,
      node_attrs: BTreeMap::new(),
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    }
  }

  #[must_use]
  pub fn with_name(mut self, name: impl Into<String>) -> Self {
    self.label = Some(NewickLabel::Name(name.into()));
    self
  }

  #[must_use]
  pub fn with_support(mut self, support: f64) -> Self {
    self.label = Some(NewickLabel::Support(support));
    self
  }

  pub fn name(&self) -> Option<&str> {
    match &self.label {
      Some(NewickLabel::Name(name)) => Some(name),
      Some(NewickLabel::Support(_)) | None => None,
    }
  }

  pub fn support(&self) -> Option<f64> {
    match &self.label {
      Some(NewickLabel::Support(support)) => Some(*support),
      Some(NewickLabel::Name(_)) | None => None,
    }
  }
}

impl Default for NewickNodeData {
  fn default() -> Self {
    Self::new()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum NewickLabel {
  Name(String),
  Support(f64),
}

impl PartialEq for NewickLabel {
  fn eq(&self, other: &Self) -> bool {
    match (self, other) {
      (Self::Name(a), Self::Name(b)) => a == b,
      (Self::Support(a), Self::Support(b)) => a.to_bits() == b.to_bits(),
      (Self::Name(_), Self::Support(_)) | (Self::Support(_), Self::Name(_)) => false,
    }
  }
}

impl Eq for NewickLabel {}

impl PartialOrd for NewickLabel {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

impl Ord for NewickLabel {
  fn cmp(&self, other: &Self) -> Ordering {
    match (self, other) {
      (Self::Name(a), Self::Name(b)) => a.cmp(b),
      (Self::Support(a), Self::Support(b)) => a.to_bits().cmp(&b.to_bits()),
      (Self::Name(_), Self::Support(_)) => Ordering::Less,
      (Self::Support(_), Self::Name(_)) => Ordering::Greater,
    }
  }
}

impl Hash for NewickLabel {
  fn hash<H: Hasher>(&self, state: &mut H) {
    std::mem::discriminant(self).hash(state);
    match self {
      Self::Name(name) => name.hash(state),
      Self::Support(support) => support.to_bits().hash(state),
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickHybrid {
  pub kind: Option<String>,
  pub index: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct NewickEdgeData {
  pub branch_length: Option<f64>,
  pub branch_attrs: BTreeMap<String, NewickValue>,
  pub raw_comments: Vec<String>,
  pub is_acceptor: bool,
}

impl NewickEdgeData {
  pub fn new() -> Self {
    Self {
      branch_length: None,
      branch_attrs: BTreeMap::new(),
      raw_comments: Vec::new(),
      is_acceptor: false,
    }
  }

  #[must_use]
  pub fn with_length(mut self, length: f64) -> Self {
    self.branch_length = Some(length);
    self
  }
}

impl Default for NewickEdgeData {
  fn default() -> Self {
    Self::new()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum NewickValue {
  Boolean(bool),
  Number(f64),
  NumberText(String),
  String(String),
  Array(Vec<NewickValue>),
}

impl PartialEq for NewickValue {
  fn eq(&self, other: &Self) -> bool {
    match (self, other) {
      (Self::Boolean(a), Self::Boolean(b)) => a == b,
      (Self::Number(a), Self::Number(b)) => a.to_bits() == b.to_bits(),
      (Self::NumberText(a), Self::NumberText(b)) | (Self::String(a), Self::String(b)) => a == b,
      (Self::Array(a), Self::Array(b)) => a == b,
      _ => false,
    }
  }
}

impl Eq for NewickValue {}

impl PartialOrd for NewickValue {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

impl Ord for NewickValue {
  fn cmp(&self, other: &Self) -> Ordering {
    match (self, other) {
      (Self::Boolean(a), Self::Boolean(b)) => a.cmp(b),
      (Self::Number(a), Self::Number(b)) => a.to_bits().cmp(&b.to_bits()),
      (Self::NumberText(a), Self::NumberText(b)) | (Self::String(a), Self::String(b)) => a.cmp(b),
      (Self::Array(a), Self::Array(b)) => a.cmp(b),
      _ => self.rank().cmp(&other.rank()),
    }
  }
}

impl NewickValue {
  fn rank(&self) -> u8 {
    match self {
      Self::Boolean(_) => 0,
      Self::Number(_) => 1,
      Self::NumberText(_) => 2,
      Self::String(_) => 3,
      Self::Array(_) => 4,
    }
  }
}

impl Hash for NewickValue {
  fn hash<H: Hasher>(&self, state: &mut H) {
    std::mem::discriminant(self).hash(state);
    match self {
      Self::Boolean(b) => b.hash(state),
      Self::Number(n) => n.to_bits().hash(state),
      Self::NumberText(s) | Self::String(s) => s.hash(state),
      Self::Array(a) => a.hash(state),
    }
  }
}
