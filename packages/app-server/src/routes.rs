use crate::args::{
  ServerAncestralArgs, ServerClockArgs, ServerMugrationArgs, ServerOptimizeArgs, ServerPruneArgs, ServerTimetreeArgs,
};
use crate::error::AppError;
use crate::sse::handle_command;
use crate::state::ServerConfig;
use app_api::datasets::discover_datasets;
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
  handle_command::<ServerAncestralArgs, _, _>(body, &config.out_dir, app_api::commands::ancestral)
}

async fn handle_clock(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerClockArgs, _, _>(body, &config.out_dir, app_api::commands::clock)
}

async fn handle_timetree(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerTimetreeArgs, _, _>(body, &config.out_dir, app_api::commands::timetree)
}

async fn handle_mugration(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerMugrationArgs, _, _>(body, &config.out_dir, app_api::commands::mugration)
}

async fn handle_optimize(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerOptimizeArgs, _, _>(body, &config.out_dir, app_api::commands::optimize)
}

async fn handle_prune(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerPruneArgs, _, _>(body, &config.out_dir, app_api::commands::prune)
}
