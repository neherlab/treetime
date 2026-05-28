use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

// --- Generic container ---

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AugurNodeDataJson<M, N> {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub generated_by: Option<AugurNodeDataJsonGeneratedBy>,

  pub nodes: BTreeMap<String, N>,

  #[serde(flatten)]
  pub metadata: M,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AugurNodeDataJsonGeneratedBy {
  pub program: String,
  pub version: String,
}

// --- Ancestral ---

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

// --- Refine (timetree) ---

pub type AugurNodeDataJsonRefine = AugurNodeDataJson<AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode>;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
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

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
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

// --- Traits (mugration) ---

pub type AugurNodeDataJsonTraits = AugurNodeDataJson<AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode>;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsMeta {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub models: Option<BTreeMap<String, AugurNodeDataJsonTraitModel>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitModel {
  pub rate: f64,
  pub alphabet: Vec<String>,
  pub equilibrium_probabilities: Vec<f64>,
  pub transition_matrix: Vec<Vec<f64>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsNode {
  #[serde(flatten)]
  pub fields: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AugurNodeDataJsonTraitsBranches {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub labels: Option<BTreeMap<String, String>>,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

// --- Annotations ---

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotations {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub nuc: Option<AugurNodeDataJsonAnnotationEntry>,

  #[serde(flatten)]
  pub other: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AugurNodeDataJsonAnnotationSegment {
  pub start: i64,
  pub end: i64,

  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
