use crate::common::AugurNodeDataJson;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type AugurNodeDataJsonRefine = AugurNodeDataJson<AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode>;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonRefineMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub alignment: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub input_tree: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub clock: Option<AugurNodeDataJsonClock>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonClock {
  pub rate: f64,
  pub intercept: f64,

  #[serde(rename = "rtt_Tmrca")]
  pub rtt_tmrca: f64,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub cov: Option<Vec<Vec<f64>>>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub rate_std: Option<f64>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonRefineNode {
  pub branch_length: f64,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<f64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub numdate: Option<f64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub clock_length: Option<f64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub mutation_length: Option<f64>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub raw_date: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub date: Option<String>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub date_inferred: Option<bool>,

  #[serde(skip_serializing_if = "Option::is_none")]
  pub num_date_confidence: Option<[f64; 2]>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}
