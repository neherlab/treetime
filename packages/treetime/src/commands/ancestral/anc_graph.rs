#![allow(clippy::default_trait_access)]
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{format_weight, nwk_read_file, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use crate::o;
use crate::seq::find_mixed_sites::MixedSite;
use eyre::Report;
use itertools::Itertools;
use maplit::{btreemap, btreeset};
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub type AncestralGraph = Graph<Node, Edge>;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Node {
  pub name: Option<String>,

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
  /// Gather all states at a given position from all child nodes
  pub fn get_letter_disambiguated(&self, pos: usize) -> BTreeSet<char> {
    self
      // Get possible states if non-consensus
      .non_consensus.get(&pos).cloned()
      // Otherwise read the sequence at that position
      .unwrap_or_else(|| btreeset!{self.seq[pos]})
  }
}

impl GraphNode for Node {}

impl Named for Node {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl NodeFromNwk for Node {
  fn from_nwk(name: Option<impl AsRef<str>>, _comments: &BTreeMap<String, String>) -> Result<Self, Report> {
    // TODO: parse mutations from comments
    Ok(Self {
      name: name.map(|n| n.as_ref().to_owned()),

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
    })
  }
}

impl NodeToNwk for Node {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations = self
      .mutations
      .iter()
      .map(|(pos, (anc, der))| format!("{anc}{}{der}", pos + 1))
      .join(",");
    BTreeMap::from([(o!("mutations"), mutations)])
  }
}

impl NodeToGraphviz for Node {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Edge {
  pub weight: Option<f64>,
}

impl GraphEdge for Edge {}

impl Weighted for Edge {
  fn weight(&self) -> Option<f64> {
    self.weight
  }
  fn set_weight(&mut self, weight: Option<f64>) {
    self.weight = weight;
  }
}

impl EdgeFromNwk for Edge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self { weight })
  }
}

impl EdgeToNwk for Edge {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight
  }
}

impl EdgeToGraphViz for Edge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight
  }
}

pub fn infer_graph() -> Result<AncestralGraph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

pub fn create_graph<P: AsRef<Path>>(tree_path: P) -> Result<AncestralGraph, Report> {
  nwk_read_file(tree_path)
}
