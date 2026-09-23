use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::fmt;
use std::hash::{DefaultHasher, Hash, Hasher};

#[derive(Clone, Debug, SmartDefault, Serialize, Deserialize)]
pub struct NewickWriteOptions {
  pub style: NwkStyle,
  pub significant_digits: Option<u8>,
  pub decimal_digits: Option<i8>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, SmartDefault, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum NwkStyle {
  Plain,
  #[default]
  Beast,
  Nhx,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NexusTree {
  pub name: String,
  pub graph: NewickGraph,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
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

  pub fn eq_ordered(&self, other: &NewickGraph) -> bool {
    if self.rooted != other.rooted {
      return false;
    }
    eq_subtree_ordered(self, self.root, other, other.root)
  }
}

impl Default for NewickGraph {
  fn default() -> Self {
    Self::new()
  }
}

impl PartialEq for NewickGraph {
  fn eq(&self, other: &Self) -> bool {
    if self.rooted != other.rooted {
      return false;
    }
    eq_subtree_unordered(self, self.root, other, other.root)
  }
}

impl Eq for NewickGraph {}

fn subtree_hash(graph: &NewickGraph, node_idx: usize) -> u64 {
  let mut hasher = DefaultHasher::new();
  let node = &graph.nodes[node_idx];
  node.name.hash(&mut hasher);
  node.confidence.map(f64::to_bits).hash(&mut hasher);
  node.node_attrs.hash(&mut hasher);
  node.raw_comments.hash(&mut hasher);
  node.hybrid.hash(&mut hasher);

  #[expect(
    clippy::collection_is_never_read,
    reason = "the vector is read by Hash::hash, which the lint does not see"
  )]
  let mut child_hashes: Vec<(u64, u64)> = node
    .children
    .iter()
    .map(|&ei| {
      let edge = &graph.edges[ei];
      let child_h = subtree_hash(graph, edge.child);
      let mut edge_hasher = DefaultHasher::new();
      edge.data.branch_length.map(f64::to_bits).hash(&mut edge_hasher);
      edge.data.branch_attrs.hash(&mut edge_hasher);
      edge.data.raw_comments.hash(&mut edge_hasher);
      edge.data.is_acceptor.hash(&mut edge_hasher);
      (child_h, edge_hasher.finish())
    })
    .collect();
  child_hashes.sort_unstable();
  child_hashes.hash(&mut hasher);
  hasher.finish()
}

fn eq_subtree_unordered(g1: &NewickGraph, n1: usize, g2: &NewickGraph, n2: usize) -> bool {
  let nd1 = &g1.nodes[n1];
  let nd2 = &g2.nodes[n2];

  if nd1.name != nd2.name
    || nd1.confidence.map(f64::to_bits) != nd2.confidence.map(f64::to_bits)
    || nd1.node_attrs != nd2.node_attrs
    || nd1.raw_comments != nd2.raw_comments
    || nd1.hybrid != nd2.hybrid
    || nd1.children.len() != nd2.children.len()
  {
    return false;
  }

  if nd1.children.is_empty() {
    return true;
  }

  let mut hashes1: Vec<(u64, u64)> = nd1
    .children
    .iter()
    .map(|&ei| {
      let edge = &g1.edges[ei];
      let child_h = subtree_hash(g1, edge.child);
      let mut eh = DefaultHasher::new();
      edge.data.branch_length.map(f64::to_bits).hash(&mut eh);
      edge.data.branch_attrs.hash(&mut eh);
      edge.data.raw_comments.hash(&mut eh);
      edge.data.is_acceptor.hash(&mut eh);
      (child_h, eh.finish())
    })
    .collect();
  let mut hashes2: Vec<(u64, u64)> = nd2
    .children
    .iter()
    .map(|&ei| {
      let edge = &g2.edges[ei];
      let child_h = subtree_hash(g2, edge.child);
      let mut eh = DefaultHasher::new();
      edge.data.branch_length.map(f64::to_bits).hash(&mut eh);
      edge.data.branch_attrs.hash(&mut eh);
      edge.data.raw_comments.hash(&mut eh);
      edge.data.is_acceptor.hash(&mut eh);
      (child_h, eh.finish())
    })
    .collect();

  hashes1.sort_unstable();
  hashes2.sort_unstable();
  hashes1 == hashes2
}

fn eq_subtree_ordered(g1: &NewickGraph, n1: usize, g2: &NewickGraph, n2: usize) -> bool {
  let nd1 = &g1.nodes[n1];
  let nd2 = &g2.nodes[n2];

  if nd1.name != nd2.name
    || nd1.confidence.map(f64::to_bits) != nd2.confidence.map(f64::to_bits)
    || nd1.node_attrs != nd2.node_attrs
    || nd1.raw_comments != nd2.raw_comments
    || nd1.hybrid != nd2.hybrid
    || nd1.children.len() != nd2.children.len()
  {
    return false;
  }

  nd1.children.iter().zip(nd2.children.iter()).all(|(&ei1, &ei2)| {
    let e1 = &g1.edges[ei1];
    let e2 = &g2.edges[ei2];
    edge_data_eq(&e1.data, &e2.data) && eq_subtree_ordered(g1, e1.child, g2, e2.child)
  })
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickEdgeEntry {
  pub parent: usize,
  pub child: usize,
  pub data: NewickEdgeData,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickNodeData {
  pub name: Option<String>,
  pub confidence: Option<f64>,
  pub node_attrs: BTreeMap<String, NewickValue>,
  pub raw_comments: Vec<String>,
  pub hybrid: Option<NewickHybrid>,
  pub children: Vec<usize>,
}

impl NewickNodeData {
  pub fn new() -> Self {
    Self {
      name: None,
      confidence: None,
      node_attrs: BTreeMap::new(),
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    }
  }

  #[must_use]
  pub fn with_name(mut self, name: impl Into<String>) -> Self {
    self.name = Some(name.into());
    self
  }
}

impl Default for NewickNodeData {
  fn default() -> Self {
    Self::new()
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NewickHybrid {
  pub kind: Option<String>,
  pub index: u32,
}

fn edge_data_eq(e1: &NewickEdgeData, e2: &NewickEdgeData) -> bool {
  let bl_eq = match (e1.branch_length, e2.branch_length) {
    (Some(a), Some(b)) => a.to_bits() == b.to_bits(),
    (None, None) => true,
    _ => false,
  };
  bl_eq && e1.branch_attrs == e2.branch_attrs && e1.raw_comments == e2.raw_comments && e1.is_acceptor == e2.is_acceptor
}

#[derive(Clone, Debug, Serialize, Deserialize)]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum NewickValue {
  Boolean(bool),
  Number(f64),
  String(String),
  Array(Vec<NewickValue>),
}

impl PartialEq for NewickValue {
  fn eq(&self, other: &Self) -> bool {
    match (self, other) {
      (Self::Boolean(a), Self::Boolean(b)) => a == b,
      (Self::Number(a), Self::Number(b)) => a.to_bits() == b.to_bits(),
      (Self::String(a), Self::String(b)) => a == b,
      (Self::Array(a), Self::Array(b)) => a == b,
      _ => false,
    }
  }
}

impl Eq for NewickValue {}

impl Hash for NewickValue {
  fn hash<H: Hasher>(&self, state: &mut H) {
    std::mem::discriminant(self).hash(state);
    match self {
      Self::Boolean(b) => b.hash(state),
      Self::Number(n) => n.to_bits().hash(state),
      Self::String(s) => s.hash(state),
      Self::Array(a) => a.hash(state),
    }
  }
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    handwritten_fmt_impl,
    reason = "array values render recursively in the Newick comment syntax"
  )
)]
impl fmt::Display for NewickValue {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      Self::Boolean(b) => write!(f, "{b}"),
      Self::Number(n) => write!(f, "{n}"),
      Self::String(s) => write!(f, "{s}"),
      Self::Array(arr) => {
        write!(f, "{{")?;
        for (i, elem) in arr.iter().enumerate() {
          if i > 0 {
            write!(f, ",")?;
          }
          write!(f, "{elem}")?;
        }
        write!(f, "}}")
      },
    }
  }
}
