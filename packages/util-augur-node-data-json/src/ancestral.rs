use crate::annotations::AugurNodeDataJsonAnnotations;
use crate::common::AugurNodeDataJson;
use deser::{Deserialize, Serialize};
use deser_value::Value;
use std::collections::BTreeMap;

pub type AugurNodeDataJsonAncestral = AugurNodeDataJson<AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode>;

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAncestralMeta {
  #[deser(skip_serializing_if = Option::is_none)]
  pub annotations: Option<AugurNodeDataJsonAnnotations>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub reference: Option<BTreeMap<String, String>>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub mask: Option<String>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAncestralNode {
  pub muts: Vec<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub sequence: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub aa_muts: Option<BTreeMap<String, Vec<String>>>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub aa_sequences: Option<BTreeMap<String, String>>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}
