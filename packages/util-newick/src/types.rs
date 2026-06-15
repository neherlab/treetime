use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::fmt;
use std::hash::{DefaultHasher, Hash, Hasher};

/// Adjacency-list representation of a Newick tree or eNewick network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickGraph {
  pub nodes: Vec<NewickNodeData>,
  pub edges: Vec<NewickEdgeEntry>,
  /// Node index of the tree root.
  pub root: usize,
  /// `[&R]` -> `Some(true)`, `[&U]` -> `Some(false)`, absent -> `None`.
  pub rooted: Option<bool>,
}

/// Edge in the adjacency list, connecting parent to child by node index.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickEdgeEntry {
  pub parent: usize,
  pub child: usize,
  pub data: NewickEdgeData,
}

/// Per-node data: name, structured annotations, raw comments, eNewick hybrid info.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickNodeData {
  pub name: Option<String>,
  /// Structured annotations from `[&...]` or `[&&NHX:...]` before `:`.
  pub node_attrs: BTreeMap<String, NewickValue>,
  /// Plain `[text]` comments without `&` prefix, preserved verbatim.
  pub raw_comments: Vec<String>,
  pub hybrid: Option<NewickHybrid>,
  /// Edge indices into `NewickGraph.edges`, in insertion (parse) order.
  pub children: Vec<usize>,
}

/// Per-edge data: branch length, structured annotations, raw comments.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NewickEdgeData {
  pub branch_length: Option<f64>,
  /// Structured annotations from after `:` (BEAST2 canonical or MrBayes position).
  pub branch_attrs: BTreeMap<String, NewickValue>,
  /// Plain `[text]` comments on the branch.
  pub raw_comments: Vec<String>,
  /// eNewick `##` (double-hash) marking this edge as the acceptor/main lineage.
  pub is_acceptor: bool,
}

/// eNewick hybrid marker parsed from `#TypeN` or `##TypeN` in node labels.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NewickHybrid {
  /// Reticulation type: `H`, `LGT`, `R`, or custom. `None` when untyped (`#1`).
  pub kind: Option<String>,
  pub index: u32,
}

/// Annotation value from BEAST `[&k=v]` or NHX `[&&NHX:k=v]` comments.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum NewickValue {
  /// BEAST `TRUE`/`FALSE` (case-insensitive).
  Boolean(bool),
  /// Bare numeric value, parsed as f64.
  Number(f64),
  /// Unquoted or quoted string value. NHX values are always String.
  String(String),
  /// BEAST `{v1,v2,v3}` array.
  Array(Vec<NewickValue>),
}

/// Output annotation style.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, SmartDefault, Serialize, Deserialize)]
pub enum NwkStyle {
  /// Name and branch length only, no annotations.
  Plain,
  /// BEAST2 `[&k=v,k2=v2]` format (default).
  #[default]
  Beast,
  /// NHX `[&&NHX:k=v:k2=v2]` format.
  Nhx,
}

/// Controls output format and float precision.
#[derive(Clone, Debug, SmartDefault, Serialize, Deserialize)]
pub struct NewickWriteOptions {
  pub style: NwkStyle,
  /// Max significant digits for branch lengths. Default: None (full precision).
  pub significant_digits: Option<u8>,
  /// Max decimal places for branch lengths. Default: None (unlimited).
  pub decimal_digits: Option<i8>,
}

/// Named tree from a Nexus TREES block.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NexusTree {
  pub name: String,
  pub graph: NewickGraph,
}

impl NewickGraph {
  /// Create an empty graph.
  pub fn new() -> Self {
    Self {
      nodes: Vec::new(),
      edges: Vec::new(),
      root: 0,
      rooted: None,
    }
  }

  /// Add a node, return its index.
  pub fn add_node(&mut self, data: NewickNodeData) -> usize {
    let idx = self.nodes.len();
    self.nodes.push(data);
    idx
  }

  /// Add an edge from parent to child, return its index.
  pub fn add_edge(&mut self, parent: usize, child: usize, data: NewickEdgeData) -> usize {
    let idx = self.edges.len();
    self.edges.push(NewickEdgeEntry { parent, child, data });
    self.nodes[parent].children.push(idx);
    idx
  }

  /// Compare two graphs with child-order sensitivity (unlike `PartialEq` which is order-insensitive).
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

impl NewickNodeData {
  /// Create a node with no name, no annotations, no children.
  pub fn new() -> Self {
    Self {
      name: None,
      node_attrs: BTreeMap::new(),
      raw_comments: Vec::new(),
      hybrid: None,
      children: Vec::new(),
    }
  }

  /// Set the node name (builder pattern).
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

impl NewickEdgeData {
  /// Create an edge with no branch length and no annotations.
  pub fn new() -> Self {
    Self {
      branch_length: None,
      branch_attrs: BTreeMap::new(),
      raw_comments: Vec::new(),
      is_acceptor: false,
    }
  }

  /// Set the branch length (builder pattern).
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

fn subtree_hash(graph: &NewickGraph, node_idx: usize) -> u64 {
  let mut hasher = DefaultHasher::new();
  let node = &graph.nodes[node_idx];
  node.name.hash(&mut hasher);
  node.node_attrs.hash(&mut hasher);
  node.raw_comments.hash(&mut hasher);
  node.hybrid.hash(&mut hasher);

  #[allow(clippy::collection_is_never_read)] // false positive: vec is read by .hash()
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

  // Compute (subtree_hash, edge_hash) for each child, sort, compare
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

fn edge_data_eq(e1: &NewickEdgeData, e2: &NewickEdgeData) -> bool {
  let bl_eq = match (e1.branch_length, e2.branch_length) {
    (Some(a), Some(b)) => a.to_bits() == b.to_bits(),
    (None, None) => true,
    _ => false,
  };
  bl_eq && e1.branch_attrs == e2.branch_attrs && e1.raw_comments == e2.raw_comments && e1.is_acceptor == e2.is_acceptor
}
