use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub struct ProgressEvent {
  pub stage: String,
  pub fraction: f64,
  pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub(crate) struct ErrorResponse {
  code: String,
  message: String,
}
