use app_api::progress::NoopProgress;
use axum::routing::post;
use axum::{Json, Router};
use serde_json::Value;

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
  tokio::task::spawn_blocking(move || app_api::commands::ancestral(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

<<<<<<< HEAD
async fn handle_clock(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || app_api::commands::clock(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_timetree(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || app_api::commands::timetree(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_mugration(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || app_api::commands::mugration(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_optimize(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || app_api::commands::optimize(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}

async fn handle_prune(Json(body): Json<Value>) -> Result<Json<Value>, (axum::http::StatusCode, String)> {
  let args = serde_json::from_value(body).map_err(json_to_axum)?;
  tokio::task::spawn_blocking(move || app_api::commands::prune(&args, &NoopProgress))
    .await
    .map_err(|e| (axum::http::StatusCode::INTERNAL_SERVER_ERROR, format!("{e}")))?
    .map_err(eyre_to_axum)?;
  Ok(Json(serde_json::json!({"status": "ok"})))
}
=======
define_handler!(
  handle_ancestral,
  app_api::TreetimeAncestralArgs,
  app_api::commands::ancestral
);
define_handler!(handle_clock, app_api::TreetimeClockArgs, app_api::commands::clock);
define_handler!(
  handle_timetree,
  app_api::TreetimeTimetreeArgs,
  app_api::commands::timetree
);
define_handler!(
  handle_mugration,
  app_api::TreetimeMugrationArgs,
  app_api::commands::mugration
);
define_handler!(
  handle_optimize,
  app_api::TreetimeOptimizeArgs,
  app_api::commands::optimize
);
define_handler!(handle_prune, app_api::TreetimePruneArgs, app_api::commands::prune);
>>>>>>> 9d9c0073 (refactor: format)
