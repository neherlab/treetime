use deser::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJson<M, N> {
  #[deser(skip_serializing_if = Option::is_none)]
  pub generated_by: Option<AugurNodeDataJsonGeneratedBy>,

  pub nodes: BTreeMap<String, N>,

  #[deser(flatten)]
  pub metadata: M,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AugurNodeDataJsonGeneratedBy {
  pub program: String,
  pub version: String,
}
