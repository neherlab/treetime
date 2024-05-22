#![allow(clippy::default_trait_access)]
use crate::graph::create_graph_from_nwk::create_graph_from_nwk_file;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named, NodeType, WithNwkComments};
use crate::o;
use crate::seq::find_mixed_sites::MixedSite;
use eyre::Report;
use itertools::Itertools;
use maplit::{btreemap, btreeset};
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::{Display, Formatter};
use std::path::Path;

pub type AncestralGraph = Graph<Node, Edge>;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,

  pub mutations: BTreeMap<usize, (char, char)>,
  pub gaps: Vec<(usize, usize)>,
  pub ambiguous: Vec<(usize, usize)>,
  pub undetermined: Vec<(usize, usize)>,
  pub mixed: Vec<MixedSite>,
  pub non_consensus: BTreeMap<usize, BTreeSet<char>>,
  pub nuc_composition: BTreeMap<char, usize>,

  pub seq: Vec<char>,

  pub subtree_profile_variable: BTreeMap<usize, Array1<f64>>,
  pub subtree_profile_fixed: BTreeMap<char, Array1<f64>>,

  pub profile_variable: BTreeMap<usize, Array1<f64>>,
  pub profile_fixed: BTreeMap<char, Array1<f64>>,

  pub outgroup_profile_variable: BTreeMap<usize, Array1<f64>>,
  pub outgroup_profile_fixed: BTreeMap<char, Array1<f64>>,

  pub expQt: Array2<f64>, // TODO: might be possible to store at the edges. Related to branch length.
}

impl Node {
  pub fn new(name: impl AsRef<str>, node_type: NodeType) -> Self {
    Self {
      name: name.as_ref().to_owned(),
      node_type,

      mutations: BTreeMap::from([]),
      gaps: vec![],
      ambiguous: vec![],
      undetermined: vec![],
      mixed: vec![],
      non_consensus: BTreeMap::new(),
      nuc_composition: BTreeMap::new(),

      seq: vec![],

      subtree_profile_variable: btreemap! {},
      subtree_profile_fixed: btreemap! {},

      profile_variable: btreemap! {},
      profile_fixed: btreemap! {},

      outgroup_profile_variable: btreemap! {},
      outgroup_profile_fixed: btreemap! {},

      expQt: Default::default(),
    }
  }

  /// Gather all states at a given position from all child nodes
  pub fn get_letter_disambiguated(&self, pos: usize) -> BTreeSet<char> {
    self
      // Get possible states if non-consensus
      .non_consensus.get(&pos).cloned()
      // Otherwise read the sequence at that position
      .unwrap_or_else(|| btreeset!{self.seq[pos]})
  }

  #[inline]
  pub const fn is_root(&self) -> bool {
    matches!(self.node_type, NodeType::Root(_))
  }

  #[inline]
  pub const fn is_internal(&self) -> bool {
    matches!(self.node_type, NodeType::Internal(_))
  }

  #[inline]
  pub const fn is_leaf(&self) -> bool {
    matches!(self.node_type, NodeType::Leaf(_))
  }
}

impl GraphNode for Node {
  fn root(name: &str) -> Self {
    Self::new(name, NodeType::Root(name.to_owned()))
  }

  fn internal(name: &str) -> Self {
    Self::new(name, NodeType::Internal(name.to_owned()))
  }

  fn leaf(name: &str) -> Self {
    Self::new(name, NodeType::Leaf(name.to_owned()))
  }

  fn set_node_type(&mut self, node_type: NodeType) {
    self.node_type = node_type;
  }
}

impl WithNwkComments for Node {
  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations = self
      .mutations
      .iter()
      .map(|(pos, (anc, der))| format!("{anc}{}{der}", pos + 1))
      .join(",");
    BTreeMap::from([(o!("mutations"), mutations)])
  }
}

impl Named for Node {
  fn name(&self) -> &str {
    &self.name
  }

  fn set_name(&mut self, name: &str) {
    self.name = name.to_owned();
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.name)
  }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Edge {
  pub weight: f64,
}

impl GraphEdge for Edge {
  fn new(weight: f64) -> Self {
    Self { weight }
  }
}

impl Weighted for Edge {
  fn weight(&self) -> f64 {
    self.weight
  }
}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

pub fn infer_graph() -> Result<AncestralGraph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

pub fn create_graph<P: AsRef<Path>>(tree_path: P) -> Result<AncestralGraph, Report> {
  create_graph_from_nwk_file::<Node, Edge>(tree_path)
}
