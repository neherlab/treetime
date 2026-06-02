use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotations {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub nuc: Option<AugurNodeDataJsonAnnotationEntry>,

  #[serde(flatten)]
  pub other: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotationEntry {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub start: Option<i64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub end: Option<i64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub strand: Option<String>,

  #[serde(rename = "type")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub entry_type: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub segments: Option<Vec<AugurNodeDataJsonAnnotationSegment>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotationSegment {
  pub start: i64,
  pub end: i64,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}
