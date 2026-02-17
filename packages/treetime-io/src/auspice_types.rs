use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::slice::Iter;
use traversal::{Bft, DftPost, DftPre};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct AuspiceGraphMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub auspice_tree_version: Option<String>,

  pub meta: AuspiceTreeMeta,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttr {
  pub value: String,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceTreeNodeAttr {
  pub fn new(value: &str) -> Self {
    Self {
      value: value.to_owned(),
      other: serde_json::Value::default(),
    }
  }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttrF64 {
  pub value: f64,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

impl AuspiceTreeNodeAttrF64 {
  pub fn new(value: f64) -> Self {
    Self {
      value,
      other: serde_json::Value::default(),
    }
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
  pub fn is_default(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Default, Serialize, Deserialize, Debug)]
pub struct AuspiceTreeNodeAttrs {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub div: Option<f64>,

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

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
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
  pub fn is_empty(&self) -> bool {
    self == &Self::default()
  }
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
pub struct StartEnd {
  pub start: isize,
  pub end: isize,

  #[serde(flatten)]
  pub other: serde_json::Value,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Segments {
  OneSegment(StartEnd),
  MultipleSegments {
    segments: Vec<StartEnd>,

    #[serde(flatten)]
    other: serde_json::Value,
  },
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
pub struct AuspiceGenomeAnnotations {
  pub nuc: AuspiceGenomeAnnotationNuc,

  #[serde(flatten)]
  pub cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,

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

#[repr(u8)]
#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum DivergenceUnits {
  NumSubstitutionsPerYearPerSite,
  #[default]
  NumSubstitutionsPerYear,
}

impl DivergenceUnits {
  ///
  /// Guesses the unit of measurement of divergence, based on the greatest value of divergence on the tree
  ///
  pub fn guess_from_max_divergence(max_divergence: f64) -> DivergenceUnits {
    // FIXME: This should be fixed upstream in augur & auspice, but it is hard to do without breaking Auspice JSON v2 format.
    // Taken from: https://github.com/nextstrain/auspice/blob/6a2d0f276fccf05bfc7084608bb0010a79086c83/src/components/tree/phyloTree/renderers.js#L376
    // A quote from there:
    //  > Prior to Jan 2020, the divergence measure was always "subs per site per year"
    //  > however certain datasets changed this to "subs per year" across entire sequence.
    //  > This distinction is not set in the JSON, so in order to correctly display the rate
    //  > we will "guess" this here. A future augur update will export this in a JSON key,
    //  > removing the need to guess
    //
    // HACK: Arbitrary threshold to make a guess
    const HACK_MAX_DIVERGENCE_THRESHOLD: f64 = 5.0;
    if max_divergence <= HACK_MAX_DIVERGENCE_THRESHOLD {
      DivergenceUnits::NumSubstitutionsPerYearPerSite
    } else {
      DivergenceUnits::NumSubstitutionsPerYear
    }
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

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AuspiceTree {
  #[serde(flatten)]
  pub data: AuspiceTreeData,

  pub tree: AuspiceTreeNode,
}

pub type AuspiceTreeNodeIter<'a> = Iter<'a, AuspiceTreeNode>;

pub type AuspiceTreeNodeIterFn<'a> = fn(&'a AuspiceTreeNode) -> AuspiceTreeNodeIter<'a>;

impl AuspiceTree {
  /// Returns iterator for breadth-first tree traversal
  pub fn iter_breadth_first<'a>(
    &'a self,
  ) -> Bft<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    Bft::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  /// Returns iterator for depth-first pre-order tree traversal
  pub fn iter_depth_first_preorder<'a>(
    &'a self,
  ) -> DftPre<'a, AuspiceTreeNode, AuspiceTreeNodeIterFn<'a>, AuspiceTreeNodeIter<'a>> {
    DftPre::new(&self.tree, |node: &'a AuspiceTreeNode| node.children.iter())
  }

  /// Returns iterator for depth-first post-order tree traversal
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

  /// Iterates over nodes and applies a function to each. Immutable version.
  pub fn map_nodes(&self, action: fn((usize, &AuspiceTreeNode))) {
    Self::map_nodes_rec(0, &self.tree, action);
  }

  fn map_nodes_mut_rec(depth: usize, node: &mut AuspiceTreeNode, action: fn((usize, &mut AuspiceTreeNode))) {
    action((depth, node));
    for child in &mut node.children {
      Self::map_nodes_mut_rec(depth + 1, child, action);
    }
  }

  /// Iterates over nodes and applies a function to each. Mutable version.
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
