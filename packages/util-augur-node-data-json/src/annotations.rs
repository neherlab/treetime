use deser::{Deserialize, Serialize};
use deser_value::Value;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotations {
  #[deser(skip_serializing_if = Option::is_none)]
  pub nuc: Option<AugurNodeDataJsonAnnotationEntry>,

  #[deser(flatten)]
  pub other: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotationEntry {
  #[deser(skip_serializing_if = Option::is_none)]
  pub start: Option<i64>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub end: Option<i64>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub strand: Option<String>,

  #[deser(rename = "type")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub entry_type: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub segments: Option<Vec<AugurNodeDataJsonAnnotationSegment>>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotationSegment {
  pub start: i64,
  pub end: i64,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}
