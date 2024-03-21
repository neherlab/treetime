use crate::alphabet::find_mutations::Mutation;
use crate::graph::create_graph_from_nwk::create_graph_from_nwk_file;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named, NodeType, WithNwkComments};
use crate::o;
use eyre::Report;
use itertools::Itertools;
use ndarray::{array, Array1, Array2};
use std::collections::BTreeMap;
use std::fmt::{Display, Formatter};
use std::path::Path;

pub type AncestralGraph = Graph<Node, Edge>;

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,
  pub joint_Lx: Array2<f64>,   // likelihood
  pub joint_Cx: Array2<usize>, // max likelihood indices
  pub mask: Option<Array1<i8>>,
  pub seq: Array1<char>,
  pub seq_ii: Array1<usize>,
  pub mutations: Vec<Mutation>,
}

impl GraphNode for Node {
  fn root(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Root(name.to_owned()),
      joint_Lx: array![[]],
      joint_Cx: array![[]],
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
      mutations: vec![],
    }
  }

  fn internal(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Internal(name.to_owned()),
      joint_Lx: array![[]],
      joint_Cx: array![[]],
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
      mutations: vec![],
    }
  }

  fn leaf(name: &str) -> Self {
    Self {
      name: name.to_owned(),
      node_type: NodeType::Leaf(name.to_owned()),
      joint_Lx: Array2::<f64>::default((0, 0)),
      joint_Cx: Array2::<usize>::default((0, 0)),
      mask: None,
      seq: Array1::<char>::default(0),
      seq_ii: Array1::<usize>::default(0),
      mutations: vec![],
    }
  }

  fn set_node_type(&mut self, node_type: NodeType) {
    self.node_type = node_type;
  }
}

impl WithNwkComments for Node {
  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations = self.mutations.iter().map(Mutation::to_string).join(",");
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

impl Node {
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

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match &self.node_type {
      NodeType::Root(weight) => write!(f, "{weight:1.4}"),
      NodeType::Internal(weight) => write!(f, "{weight:1.4}"),
      NodeType::Leaf(name) => write!(f, "{name}"),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
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
