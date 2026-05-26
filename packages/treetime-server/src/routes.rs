use axum::routing::post;
use axum::{Json, Router};
use serde_json::Value;
use treetime_api::progress::NoopProgress;

fn eyre_to_axum(err: eyre::Report) -> (axum::http::StatusCode, String) {
  (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{err:?}"))
}

fn json_to_axum(err: serde_json::Error) -> (axum::http::StatusCode, String) {
  (axum::http::StatusCode::BAD_REQUEST, format!("{err}"))
}

pub fn api_routes() -> Router {
  Router::new()
    .route("/ancestral", post(handle_ancestral))
    .route("/clock", post(handle_clock))
    .route("/timetree", post(handle_timetree))
    .route("/mugration", post(handle_mugration))
    .route("/optimize", post(handle_optimize))
    .route("/prune", post(handle_prune))
}

async fn handle_ancestral(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || treetime_api::commands::ancestral(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_clock(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || treetime_api::commands::clock(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_timetree(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || treetime_api::commands::timetree(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_mugration(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || treetime_api::commands::mugration(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_optimize(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || treetime_api::commands::optimize(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_prune(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || treetime_api::commands::prune(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}
