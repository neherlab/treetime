#![allow(
  clippy::empty_structs_with_brackets,
  reason = "ToSchema types serialize as an empty JSON object in the OpenAPI contract; a unit struct would change the generated schema"
)]

use serde_json::Value;
use utoipa::ToSchema;

#[derive(ToSchema)]
pub struct VersionInfo {
  pub version: String,
}

#[derive(ToSchema)]
pub struct DatasetInfo {
  pub name: String,
  pub files: Vec<String>,
}

#[derive(ToSchema)]
pub struct ProgressEvent {
  pub stage: String,
  pub fraction: f64,
  pub message: String,
}

#[derive(ToSchema)]
pub enum LogLevel {
  Trace,
  Debug,
  Info,
  Warn,
  Error,
}

#[derive(ToSchema)]
pub struct LogEvent {
  pub level: LogLevel,
  pub message: String,
}

#[derive(ToSchema)]
pub struct ErrorResponse {
  pub code: String,
  pub message: String,
}

#[derive(ToSchema)]
pub struct AncestralResult {
  pub model_name: String,
}

#[derive(ToSchema)]
pub struct ClockResult {
  #[schema(value_type = Object)]
  pub clock_model: Value,
  #[schema(value_type = Vec<Object>)]
  pub regression_results: Vec<Value>,
}

#[derive(ToSchema)]
pub struct TimetreeResult {
  #[schema(value_type = Object)]
  pub clock_model: Value,
}

#[derive(ToSchema)]
pub struct MugrationResult {
  pub log_lh: f64,
}

#[derive(ToSchema)]
pub struct OptimizeResult {}

#[derive(ToSchema)]
pub struct PruneResult {}
