use crate::common::AugurNodeDataJson;
use deser::{Deserialize, Serialize};
use deser_value::Value;
use std::collections::BTreeMap;

pub type AugurNodeDataJsonRefine = AugurNodeDataJson<AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode>;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonRefineMeta {
  #[deser(skip_serializing_if = Option::is_none)]
  pub alignment: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub input_tree: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub clock: Option<AugurNodeDataJsonClock>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonClock {
  pub rate: f64,
  pub intercept: f64,

  #[deser(rename = "rtt_Tmrca")]
  pub rtt_tmrca: f64,

  #[deser(skip_serializing_if = Option::is_none)]
  pub cov: Option<Vec<Vec<f64>>>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub rate_std: Option<f64>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonRefineNode {
  pub branch_length: f64,

  #[deser(skip_serializing_if = Option::is_none)]
  pub numdate: Option<f64>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub clock_length: Option<f64>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub mutation_length: Option<f64>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub raw_date: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub date: Option<String>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub date_inferred: Option<bool>,

  #[deser(skip_serializing_if = Option::is_none)]
  pub num_date_confidence: Option<[f64; 2]>,

  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}
