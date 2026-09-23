use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::slice::Iter;
use traversal::{Bft, DftPost, DftPre};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct AuspiceGraphMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  auspice_tree_version: Option<String>,

  meta: AuspiceTreeMeta,

  #[serde(flatten)]
  other: serde_json::Value,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttrF64 {
  value: f64,

  #[serde(flatten)]
  other: serde_json::Value,
}

impl AuspiceTreeNodeAttrF64 {
  pub fn new(value: f64) -> Self {
    Self {
      value,
      other: serde_json::Value::default(),
    }
  }
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum DivergenceUnits {
  NumSubstitutionsPerYearPerSite,
  #[default]
  NumSubstitutionsPerYear,
}

impl DivergenceUnits {
  pub fn guess_from_max_divergence(max_divergence: f64) -> DivergenceUnits {
    const HACK_MAX_DIVERGENCE_THRESHOLD: f64 = 5.0;
    if max_divergence <= HACK_MAX_DIVERGENCE_THRESHOLD {
      DivergenceUnits::NumSubstitutionsPerYearPerSite
    } else {
      DivergenceUnits::NumSubstitutionsPerYear
    }
  }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTree {
  #[serde(flatten)]
  pub data: AuspiceTreeData,

  pub tree: AuspiceTreeNode,
}

impl AuspiceTree {
  pub fn iter_breadth_first<'a>(
    &'a self,
  ) -> Bft<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    Bft::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  pub fn iter_depth_first_preorder<'a>(
    &'a self,
  ) -> DftPre<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    DftPre::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  pub fn iter_depth_first_postorder<'a>(
    &'a self,
  ) -> DftPost<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    DftPost::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  fn map_nodes_rec(depth: usize, node: &AuspiceTreeNode, action: fn((usize, &AuspiceTreeNode))) {
    action((depth, node));
    for child in &node.children {
      Self::map_nodes_rec(depth + 1, child, action);
    }
  }

  pub fn map_nodes(&self, action: fn((usize, &AuspiceTreeNode))) {
    Self::map_nodes_rec(0, &self.tree, action);
  }

  fn map_nodes_mut_rec(depth: usize, node: &mut AuspiceTreeNode, action: fn((usize, &mut AuspiceTreeNode))) {
    action((depth, node));
    for child in &mut node.children {
      Self::map_nodes_mut_rec(depth + 1, child, action);
    }
  }

  pub fn map_nodes_mut(&mut self, action: fn((usize, &mut AuspiceTreeNode))) {
    Self::map_nodes_mut_rec(0, &mut self.tree, action);
  }

  pub fn root_sequence(&self) -> Option<&str> {
    self
      .data
      .root_sequence
      .as_ref()
      .and_then(|root_sequence| root_sequence.get("nuc"))
      .map(String::as_str)
  }
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeData {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub version: Option<String>,

  pub meta: AuspiceTreeMeta,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub root_sequence: Option<BTreeMap<String, String>>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeMeta {
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub title: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub description: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub updated: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub genome_annotations: Option<AuspiceGenomeAnnotations>,

  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub colorings: Vec<AuspiceColoring>,

  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub panels: Vec<String>,

  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub filters: Vec<String>,

  #[serde(default, skip_serializing_if = "AuspiceDisplayDefaults::is_empty")]
  pub display_defaults: AuspiceDisplayDefaults,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub geo_resolutions: Option<serde_json::Value>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuspiceColoring {
  #[serde(rename = "type")]
  pub type_: String,

  pub key: String,

  pub title: String,

  #[serde(skip_serializing_if = "Vec::is_empty")]
  #[serde(default)]
  pub scale: Vec<[String; 2]>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Serialize, Deserialize, Eq, PartialEq, Debug)]
pub struct AuspiceDisplayDefaults {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub branch_label: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub color_by: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub distance_measure: Option<String>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceDisplayDefaults {
  fn is_empty(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotations {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub nuc: Option<AuspiceGenomeAnnotationNuc>,

  #[serde(flatten)]
  pub cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotationNuc {
  pub start: isize,

  pub end: isize,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub strand: Option<String>,

  #[serde(rename = "type", default, skip_serializing_if = "Option::is_none")]
  pub r#type: Option<String>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotationCds {
  #[serde(rename = "type", default, skip_serializing_if = "Option::is_none")]
  pub r#type: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub gene: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub color: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub display_name: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub description: Option<String>,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub strand: Option<String>,

  #[serde(flatten)]
  pub segments: Segments,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(untagged)]
#[serde(rename_all = "kebab-case")]
pub enum Segments {
  OneSegment(StartEnd),
  MultipleSegments {
    segments: Vec<StartEnd>,

    #[serde(flatten)]
    other: serde_json::Value,
  },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StartEnd {
  pub start: isize,
  pub end: isize,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

pub type AuspiceTreeNodeIterFn<'a> = fn(&'a AuspiceTreeNode) -> AuspiceTreeNodeIter<'a>;

pub type AuspiceTreeNodeIter<'a> = Iter<'a, AuspiceTreeNode>;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNode {
  pub name: String,

  #[serde(default, skip_serializing_if = "AuspiceTreeBranchAttrs::is_default")]
  pub branch_attrs: AuspiceTreeBranchAttrs,

  pub node_attrs: AuspiceTreeNodeAttrs,

  #[serde(skip_serializing_if = "Vec::is_empty")]
  #[serde(default)]
  pub children: Vec<AuspiceTreeNode>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Eq, PartialEq, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeBranchAttrs {
  pub mutations: BTreeMap<String, Vec<String>>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub labels: Option<AuspiceTreeBranchAttrsLabels>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceTreeBranchAttrs {
  #[inline]
  fn is_default(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Default, Eq, PartialEq, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeBranchAttrsLabels {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub aa: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub clade: Option<String>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttrs {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub div: Option<f64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub num_date: Option<AuspiceNumDate>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub bad_branch: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub clade_membership: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub region: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub country: Option<AuspiceTreeNodeAttr>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub division: Option<AuspiceTreeNodeAttr>,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttr {
  value: String,

  #[serde(flatten)]
  other: serde_json::Value,
}

impl AuspiceTreeNodeAttr {
  pub fn new(value: &str) -> Self {
    Self {
      value: value.to_owned(),
      other: serde_json::Value::default(),
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceNumDate {
  pub value: f64,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<[f64; 2]>,
}
