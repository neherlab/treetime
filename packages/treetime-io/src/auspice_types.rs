use deser::{Deserialize, Serialize};
use deser_value::{Map, Value, from_value};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;

#[repr(u8)]
#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub enum DivergenceUnits {
  NumSubstitutionsPerYearPerSite,
  #[default]
  NumSubstitutionsPerYear,
}

impl DivergenceUnits {}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceTree {
  #[deser(flatten)]
  pub data: AuspiceTreeData,

  pub tree: AuspiceTreeNode,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeData {
  #[deser(skip_serializing_if = Option::is_none)]
  pub version: Option<String>,

  pub meta: AuspiceTreeMeta,

  #[deser(skip_serializing_if = Option::is_none)]
  pub root_sequence: Option<BTreeMap<String, String>>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeMeta {
  #[deser(default, skip_serializing_if = Option::is_none)]
  pub title: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub description: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub updated: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub genome_annotations: Option<AuspiceGenomeAnnotations>,

  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub colorings: Vec<AuspiceColoring>,

  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub panels: Vec<String>,

  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub filters: Vec<String>,

  #[deser(default, skip_serializing_if = AuspiceDisplayDefaults::is_empty)]
  pub display_defaults: AuspiceDisplayDefaults,

  #[deser(skip_serializing_if = Option::is_none)]
  pub geo_resolutions: Option<Value>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuspiceColoring {
  #[deser(rename = "type")]
  pub type_: String,

  pub key: String,

  pub title: String,

  #[deser(skip_serializing_if = Vec::is_empty)]
  #[deser(default)]
  pub scale: Vec<[String; 2]>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Default, Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct AuspiceDisplayDefaults {
  #[deser(skip_serializing_if = Option::is_none)]
  pub branch_label: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub color_by: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub distance_measure: Option<String>,

  #[deser(flatten)]
  pub other: Map,
}

impl AuspiceDisplayDefaults {
  fn is_empty(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotations {
  #[deser(skip_serializing_if = Option::is_none)]
  pub nuc: Option<AuspiceGenomeAnnotationNuc>,

  #[deser(flatten)]
  pub cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotationNuc {
  pub start: isize,

  pub end: isize,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub strand: Option<String>,

  #[deser(rename = "type", default, skip_serializing_if = Option::is_none)]
  pub r#type: Option<String>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceGenomeAnnotationCds {
  #[deser(rename = "type", default, skip_serializing_if = Option::is_none)]
  pub r#type: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub gene: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub color: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub display_name: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub description: Option<String>,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub strand: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub start: Option<isize>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub end: Option<isize>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub segments: Option<Vec<StartEnd>>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StartEnd {
  pub start: isize,
  pub end: isize,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeNode {
  pub name: String,

  #[deser(default, skip_serializing_if = AuspiceTreeBranchAttrs::is_default)]
  pub branch_attrs: AuspiceTreeBranchAttrs,

  pub node_attrs: AuspiceTreeNodeAttrs,

  #[deser(skip_serializing_if = Vec::is_empty)]
  #[deser(default)]
  pub children: Vec<AuspiceTreeNode>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Default, Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeBranchAttrs {
  pub mutations: BTreeMap<String, Vec<String>>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub labels: Option<AuspiceTreeBranchAttrsLabels>,

  #[deser(flatten)]
  pub other: Map,
}

impl AuspiceTreeBranchAttrs {
  #[inline]
  fn is_default(&self) -> bool {
    self == &Self::default()
  }
}

#[derive(Clone, Default, Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeBranchAttrsLabels {
  #[deser(skip_serializing_if = Option::is_none)]
  pub aa: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub clade: Option<String>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeNodeAttrs {
  #[deser(skip_serializing_if = Option::is_none)]
  pub div: Option<f64>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub num_date: Option<AuspiceNumDate>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub bad_branch: Option<AuspiceTreeNodeAttr>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub clade_membership: Option<AuspiceTreeNodeAttr>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub region: Option<AuspiceTreeNodeAttr>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub country: Option<AuspiceTreeNodeAttr>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub division: Option<AuspiceTreeNodeAttr>,

  #[deser(flatten)]
  pub other: Map,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceTreeNodeAttr {
  value: String,

  #[deser(flatten)]
  other: Map,
}

impl AuspiceTreeNodeAttr {
  pub fn new(value: &str) -> Self {
    Self {
      value: value.to_owned(),
      other: Map::new(),
    }
  }

  pub fn value(&self) -> &str {
    &self.value
  }

  pub fn state_confidence(&self) -> Result<Option<BTreeMap<String, f64>>, Report> {
    match self.other.get("confidence") {
      Some(confidence) if confidence.is_map() => from_value(confidence)
        .map(Some)
        .wrap_err("When reading the state probabilities of a node attribute"),
      _ => Ok(None),
    }
  }
}

impl AuspiceTreeNodeAttrs {
  pub fn attr(&self, key: &str) -> Result<Option<AuspiceTreeNodeAttr>, Report> {
    let named = match key {
      "bad_branch" => &self.bad_branch,
      "clade_membership" => &self.clade_membership,
      "region" => &self.region,
      "country" => &self.country,
      "division" => &self.division,
      _ => {
        return self
          .other
          .get(key)
          .map(from_value)
          .transpose()
          .wrap_err_with(|| format!("When reading the node attribute `{key}`"));
      },
    };
    Ok(named.clone())
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AuspiceNumDate {
  pub value: f64,

  #[deser(skip_serializing_if = Option::is_none)]
  pub confidence: Option<[f64; 2]>,
}
