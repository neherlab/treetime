use crate::common::AugurNodeDataJson;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type AugurNodeDataJsonTraits = AugurNodeDataJson<AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode>;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub models: Option<BTreeMap<String, AugurNodeDataJsonTraitModel>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitModel {
  pub rate: f64,
  pub alphabet: Vec<String>,
  pub equilibrium_probabilities: Vec<f64>,
  pub transition_matrix: Vec<Vec<f64>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsNode {
  #[serde(flatten)]
  pub fields: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsBranches {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub labels: Option<BTreeMap<String, String>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}
