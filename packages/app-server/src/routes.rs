use crate::commands::ancestral::{AncestralArgs, run_ancestral};
use crate::commands::clock::{ClockArgs, run_clock};
use crate::commands::mugration::{MugrationArgs, run_mugration};
use crate::commands::optimize::{OptimizeArgs, run_optimize};
use crate::commands::prune::{PruneArgs, run_prune};
use crate::commands::timetree::{TimetreeArgs, run_timetree};
use crate::error::AppError;
use crate::sse::handle_command;
use crate::state::ServerConfig;
use app_datasets::discover_datasets;
use axum::extract::State;
use axum::response::Response;
use axum::routing::{get, post};
use axum::{Json, Router};
use serde_json::Value;
use std::sync::Arc;
use treetime_schema::version_info;

pub fn api_routes(config: ServerConfig) -> Router {
  let state = Arc::new(config);
  Router::new()
    .route("/health", get(handle_health))
    .route("/version", get(handle_version))
    .route("/datasets", get(handle_datasets))
    .route("/ancestral", post(handle_ancestral))
    .route("/clock", post(handle_clock))
    .route("/timetree", post(handle_timetree))
    .route("/mugration", post(handle_mugration))
    .route("/optimize", post(handle_optimize))
    .route("/prune", post(handle_prune))
    .with_state(state)
}

async fn handle_health() -> Json<Value> {
  Json(serde_json::json!({
    "status": "ok",
    "version": version_info().version,
  }))
}

async fn handle_version() -> Result<Json<Value>, AppError> {
  let value = serde_json::to_value(version_info())?;
  Ok(Json(value))
}

async fn handle_datasets(State(config): State<Arc<ServerConfig>>) -> Result<Json<Value>, AppError> {
  let datasets = discover_datasets(&config.data_dir);
  let value = serde_json::to_value(datasets)?;
  Ok(Json(value))
}

async fn handle_ancestral(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<AncestralArgs, _>(body, &config.out_dir, run_ancestral)
}

async fn handle_clock(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ClockArgs, _>(body, &config.out_dir, run_clock)
}

async fn handle_timetree(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<TimetreeArgs, _>(body, &config.out_dir, run_timetree)
}

async fn handle_mugration(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<MugrationArgs, _>(body, &config.out_dir, run_mugration)
}

async fn handle_optimize(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<OptimizeArgs, _>(body, &config.out_dir, run_optimize)
}

async fn handle_prune(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<PruneArgs, _>(body, &config.out_dir, run_prune)
}
