use crate::commands::ancestral::{AncestralArgs, run_ancestral};
use crate::commands::clock::{ClockArgs, run_clock};
use crate::commands::mugration::{MugrationArgs, run_mugration};
use crate::commands::optimize::{OptimizeArgs, run_optimize};
use crate::commands::prune::{PruneArgs, run_prune};
use crate::commands::timetree::{TimetreeArgs, run_timetree};
use crate::contract::{
  AncestralResult, ClockResult, DatasetInfo, ErrorResponse, LogEvent, MugrationResult, OptimizeResult, ProgressEvent,
  PruneResult, TimetreeResult, VersionInfo,
};
use crate::error::AppError;
use crate::sse::handle_command;
use crate::state::ServerConfig;
use app_datasets::discover_datasets;
use axum::extract::State;
use axum::response::Response;
use axum::routing::get;
use axum::{Json, Router};
use serde_json::Value;
use std::sync::Arc;
use treetime_schema::version_info;
use utoipa::OpenApi;
use utoipa::openapi::extensions::Extensions;
use utoipa::openapi::{ContactBuilder, LicenseBuilder};
use utoipa_axum::router::OpenApiRouter;
use utoipa_axum::routes;

pub fn api_routes(config: ServerConfig) -> Router {
  let (router, _api) = api_router().with_state(Arc::new(config)).split_for_parts();
  router
    .route("/api/health", get(handle_health))
    .route("/api/openapi.json", get(handle_openapi))
}

pub fn api_doc() -> utoipa::openapi::OpenApi {
  let mut api = api_router().to_openapi();
  api.info.title = "TreeTime API".to_owned();
  api.info.version = env!("CARGO_PKG_VERSION").to_owned();
  api.info.description = Some(env!("CARGO_PKG_DESCRIPTION").to_owned());
  api.info.contact = Some(
    ContactBuilder::new()
      .name(Some("NeherLab"))
      .url(Some(env!("CARGO_PKG_HOMEPAGE")))
      .build(),
  );
  api.info.license = Some(LicenseBuilder::new().name(env!("CARGO_PKG_LICENSE")).build());
  api.merge(SharedSchemas::openapi());
  set_bridge_types(&mut api);
  api
}

fn api_router() -> OpenApiRouter<Arc<ServerConfig>> {
  OpenApiRouter::new()
    .routes(routes!(handle_version))
    .routes(routes!(handle_datasets))
    .routes(routes!(handle_ancestral))
    .routes(routes!(handle_clock))
    .routes(routes!(handle_timetree))
    .routes(routes!(handle_mugration))
    .routes(routes!(handle_optimize))
    .routes(routes!(handle_prune))
}

fn set_bridge_types(api: &mut utoipa::openapi::OpenApi) {
  for item in api.paths.paths.values_mut() {
    if let Some(op) = item.get.as_mut() {
      op.extensions = Some(Extensions::from_iter([(
        "x-bridge-type",
        Value::String("query".to_owned()),
      )]));
    }
    if let Some(op) = item.post.as_mut() {
      op.extensions = Some(Extensions::from_iter([(
        "x-bridge-type",
        Value::String("command".to_owned()),
      )]));
    }
  }
}

async fn handle_health() -> Json<Value> {
  Json(serde_json::json!({
    "status": "ok",
    "version": version_info().version,
  }))
}

async fn handle_openapi() -> Result<Json<Value>, AppError> {
  Ok(Json(serde_json::to_value(api_doc())?))
}

#[utoipa::path(
  get,
  path = "/api/version",
  operation_id = "version",
  responses((status = 200, description = "Version information", body = VersionInfo))
)]
async fn handle_version() -> Result<Json<Value>, AppError> {
  let value = serde_json::to_value(version_info())?;
  Ok(Json(value))
}

#[utoipa::path(
  get,
  path = "/api/datasets",
  operation_id = "datasets",
  responses((status = 200, description = "Available datasets", body = Vec<DatasetInfo>))
)]
async fn handle_datasets(State(config): State<Arc<ServerConfig>>) -> Result<Json<Value>, AppError> {
  let datasets = discover_datasets(&config.data_dir);
  let value = serde_json::to_value(datasets)?;
  Ok(Json(value))
}

#[utoipa::path(
  post,
  path = "/api/ancestral",
  operation_id = "ancestral",
  request_body = AncestralArgs,
  responses((status = 200, description = "SSE stream with progress events and final result", body = AncestralResult, content_type = "text/event-stream"))
)]
async fn handle_ancestral(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<AncestralArgs, _>(body, &config.out_dir, run_ancestral)
}

#[utoipa::path(
  post,
  path = "/api/clock",
  operation_id = "clock",
  request_body = ClockArgs,
  responses((status = 200, description = "SSE stream with progress events and final result", body = ClockResult, content_type = "text/event-stream"))
)]
async fn handle_clock(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ClockArgs, _>(body, &config.out_dir, run_clock)
}

#[utoipa::path(
  post,
  path = "/api/timetree",
  operation_id = "timetree",
  request_body = TimetreeArgs,
  responses((status = 200, description = "SSE stream with progress events and final result", body = TimetreeResult, content_type = "text/event-stream"))
)]
async fn handle_timetree(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<TimetreeArgs, _>(body, &config.out_dir, run_timetree)
}

#[utoipa::path(
  post,
  path = "/api/mugration",
  operation_id = "mugration",
  request_body = MugrationArgs,
  responses((status = 200, description = "SSE stream with progress events and final result", body = MugrationResult, content_type = "text/event-stream"))
)]
async fn handle_mugration(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<MugrationArgs, _>(body, &config.out_dir, run_mugration)
}

#[utoipa::path(
  post,
  path = "/api/optimize",
  operation_id = "optimize",
  request_body = OptimizeArgs,
  responses((status = 200, description = "SSE stream with progress events and final result", body = OptimizeResult, content_type = "text/event-stream"))
)]
async fn handle_optimize(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<OptimizeArgs, _>(body, &config.out_dir, run_optimize)
}

#[utoipa::path(
  post,
  path = "/api/prune",
  operation_id = "prune",
  request_body = PruneArgs,
  responses((status = 200, description = "SSE stream with progress events and final result", body = PruneResult, content_type = "text/event-stream"))
)]
async fn handle_prune(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<PruneArgs, _>(body, &config.out_dir, run_prune)
}

#[derive(OpenApi)]
#[openapi(components(schemas(ProgressEvent, LogEvent, ErrorResponse)))]
struct SharedSchemas;
