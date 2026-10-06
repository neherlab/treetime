use deser::{Deserialize, Serialize};
use schemars::JsonSchema;

#[derive(Debug, Clone, JsonSchema, Serialize, Deserialize)]
pub struct ProgressEvent {
  pub stage: String,
  pub fraction: f64,
  pub message: String,
}

#[derive(Debug, Clone, JsonSchema, Serialize, Deserialize)]
pub(crate) struct ErrorResponse {
  code: String,
  message: String,
}
