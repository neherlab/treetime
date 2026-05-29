use crate::annotations::AugurNodeDataJsonAnnotations;
use crate::common::AugurNodeDataJson;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type AugurNodeDataJsonAncestral = AugurNodeDataJson<AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode>;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAncestralMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub annotations: Option<AugurNodeDataJsonAnnotations>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub reference: Option<BTreeMap<String, String>>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub mask: Option<String>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAncestralNode {
  pub muts: Vec<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub sequence: Option<String>,

  // AA reconstruction (future: kb/proposals/node-data-json-aa-reconstruction.md)
  #[serde(skip_serializing_if = "Option::is_none")]
  pub aa_muts: Option<BTreeMap<String, Vec<String>>>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub aa_sequences: Option<BTreeMap<String, String>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}
