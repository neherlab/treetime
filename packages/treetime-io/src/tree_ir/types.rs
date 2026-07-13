use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use crate::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};
use crate::tree_ir::mutation::{TreeIrIndel, TreeIrSub};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::node::{GraphNode, Named};

/// Format-neutral tree node.
///
/// Carries the union of node-level data the supported tree formats can represent.
/// Each adapter projects the subset its format supports; data a format cannot carry
/// is dropped (with a warning where it is lossy and silent omission would surprise).
///
/// `div` is cumulative divergence from the root (root = 0), precomputed at projection
/// time so every adapter emits a consistent value without re-deriving it. The
/// projection chooses the per-branch contribution source (mutation counts when
/// divergence is measured in mutations, branch length otherwise), matching augur's
/// `node_div` preference order.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TreeIrNode {
  pub name: Option<String>,
  pub desc: Option<String>,
  pub div: Option<f64>,
  /// Numeric (decimal-year) date or inferred time.
  pub date: Option<f64>,
  /// Date confidence interval `[lower, upper]`.
  pub date_confidence: Option<[f64; 2]>,
  /// Whether the date was inferred (internal nodes) rather than observed (tips).
  pub date_is_inferred: bool,
  /// Original date string for tips (augur `num_date.raw_value`).
  pub date_raw_value: Option<String>,
  /// Whether the node is flagged as a bad branch / temporal outlier.
  pub bad_branch: bool,
  /// Discrete trait assignments keyed by attribute name (e.g. mugration `country`).
  pub traits: BTreeMap<String, TreeIrTrait>,
}

impl GraphNode for TreeIrNode {}

impl Named for TreeIrNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|name| name.as_ref().to_owned());
  }
}

impl NodeToNwk for TreeIrNode {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

impl NodeFromNwk for TreeIrNode {
  fn from_nwk(
    name: Option<impl AsRef<str>>,
    _confidence: Option<f64>,
    _comments: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|name| name.as_ref().to_owned()),
      ..Self::default()
    })
  }
}

impl NodeToGraphviz for TreeIrNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

/// Discrete trait assignment for a node.
///
/// `value` is the assigned state; `confidence` maps candidate states to their
/// marginal probabilities; `entropy` is the Shannon entropy of that distribution.
/// Mirrors augur's `node_attrs.<trait>` object (`value`, `confidence`, `entropy`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TreeIrTrait {
  pub value: String,
  #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
  pub confidence: BTreeMap<String, f64>,
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub entropy: Option<f64>,
}

/// Format-neutral tree edge (branch leading to the child node).
///
/// `mutations` and `indels` describe the events on the branch above the child.
/// `branch_length` is in substitutions per site. `gamma` is the branch-rate
/// multiplier from timetree relaxed-clock inference.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TreeIrEdge {
  pub branch_length: Option<f64>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub mutations: Vec<TreeIrSub>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub indels: Vec<TreeIrIndel>,
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub gamma: Option<f64>,
}

impl GraphEdge for TreeIrEdge {}

impl HasBranchLength for TreeIrEdge {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, branch_length: Option<f64>) {
    self.branch_length = branch_length;
  }
}

impl EdgeToNwk for TreeIrEdge {
  fn nwk_weight(&self) -> Option<f64> {
    self.branch_length
  }
}

impl EdgeFromNwk for TreeIrEdge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length: weight,
      ..Self::default()
    })
  }
}

impl EdgeToGraphviz for TreeIrEdge {}

/// Format-neutral graph-level data.
///
/// Carries metadata that does not belong to any single node: dataset title and
/// description, per-gene root reference sequences, the set of discrete trait
/// attributes present (to build colorings), and flags recording which node/branch
/// data classes are populated (to decide which colorings and properties to emit).
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct TreeIrData {
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub title: Option<String>,
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub description: Option<String>,
  /// Root reference sequence per gene track (e.g. `"nuc" -> "ACGT..."`).
  #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
  pub root_sequence: BTreeMap<String, String>,
  /// Names of discrete trait attributes present on nodes.
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub trait_attrs: Vec<String>,
  /// Whether nodes carry dates (drives the `num_date` coloring and node attribute).
  pub has_dates: bool,
  /// Whether nodes carry the bad-branch flag (drives the `bad_branch` coloring).
  pub has_bad_branch: bool,
  /// Whether branches carry mutations (drives the genotype `gt` coloring).
  pub has_mutations: bool,
}
