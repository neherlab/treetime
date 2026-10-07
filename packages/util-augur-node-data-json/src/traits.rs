use crate::common::AugurNodeDataJson;
use deser::{Deserialize, Serialize};
use deser_value::Value;
use std::collections::BTreeMap;

pub type AugurNodeDataJsonTraits = AugurNodeDataJson<AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode>;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsMeta {
  #[deser(skip_serializing_if = Option::is_none)]
  pub models: Option<BTreeMap<String, AugurNodeDataJsonTraitModel>>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitModel {
  pub rate: f64,
  pub alphabet: Vec<String>,
  pub equilibrium_probabilities: Vec<f64>,
  pub transition_matrix: Vec<Vec<f64>>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsNode {
  #[deser(flatten)]
  pub fields: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct AugurNodeDataJsonTraitsBranches {
  #[deser(skip_serializing_if = Option::is_none)]
  pub labels: Option<BTreeMap<String, String>>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}
