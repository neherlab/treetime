<<<<<<< HEAD
use app_api::progress::NoopProgress;
<<<<<<< HEAD
=======
use crate::args::{
  ServerAncestralArgs, ServerClockArgs, ServerMugrationArgs, ServerOptimizeArgs, ServerPruneArgs, ServerTimetreeArgs,
};
use crate::error::AppError;
use crate::sse::handle_command;
use crate::state::ServerConfig;
use app_api::datasets::discover_datasets;
use app_api::version::version_info;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)
=======
use app_api::version::version_info;
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
use axum::http::StatusCode;
use axum::routing::{get, post};
use axum::{Json, Router};
<<<<<<< HEAD
use serde_json::Value;

fn error_response(status: StatusCode, code: &str, message: String) -> (StatusCode, Json<Value>) {
  (status, Json(serde_json::json!({"code": code, "message": message})))
}

fn eyre_to_response(err: eyre::Report) -> (StatusCode, Json<Value>) {
  error_response(
    StatusCode::INTERNAL_SERVER_ERROR,
    "computation_failure",
    format!("{err:?}"),
  )
}

fn json_to_response(err: serde_json::Error) -> (StatusCode, Json<Value>) {
  error_response(StatusCode::BAD_REQUEST, "invalid_input", format!("{err}"))
}

fn join_to_response(err: tokio::task::JoinError) -> (StatusCode, Json<Value>) {
  error_response(StatusCode::INTERNAL_SERVER_ERROR, "internal_error", format!("{err}"))
}
=======
use eyre::Report;
use serde::de::DeserializeOwned;
use serde_json::Value;
use std::convert::Infallible;
use tokio::sync::mpsc;
use tokio_stream::StreamExt as _;
use tokio_stream::wrappers::UnboundedReceiverStream;
use treetime::progress::ProgressEvent;
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)
=======
=======
use axum::extract::State;
>>>>>>> 6fc31936 (feat(server): add CLI args for jobs, data dir, and output dir)
use axum::response::Response;
use axum::routing::{get, post};
use axum::{Json, Router};
use serde_json::Value;
<<<<<<< HEAD
>>>>>>> 759dbc26 (refactor(app-server): extract error, consolidate SSE transport into modules)
=======
use std::sync::Arc;
>>>>>>> 6fc31936 (feat(server): add CLI args for jobs, data dir, and output dir)

<<<<<<< HEAD
<<<<<<< HEAD
pub fn api_routes(config: ServerConfig) -> Router {
  let state = Arc::new(config);
=======
async fn handle_version() -> Json<Value> {
  Json(serde_json::to_value(version_info()).unwrap_or_default())
=======
async fn handle_version() -> Result<Json<Value>, (StatusCode, Json<Value>)> {
<<<<<<< HEAD
  serde_json::to_value(version_info())
    .map(Json)
    .map_err(|err| error_response(StatusCode::INTERNAL_SERVER_ERROR, "serialization_error", format!("{err}")))
>>>>>>> b18af097 (fix(app-server): use structured error format and fix unwrap_or_default)
=======
  serde_json::to_value(version_info()).map(Json).map_err(|err| {
    error_response(
      StatusCode::INTERNAL_SERVER_ERROR,
      "serialization_error",
      format!("{err}"),
    )
  })
>>>>>>> f1127239 (refactor: format)
}

pub fn api_routes() -> Router {
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
  Router::new()
    .route("/version", get(handle_version))
    .route("/ancestral", post(handle_ancestral))
    .route("/clock", post(handle_clock))
    .route("/timetree", post(handle_timetree))
    .route("/mugration", post(handle_mugration))
    .route("/optimize", post(handle_optimize))
    .route("/prune", post(handle_prune))
    .with_state(state)
}

<<<<<<< HEAD
<<<<<<< HEAD
macro_rules! define_handler {
  ($handler_name:ident, $args_type:ty, $api_fn:path) => {
    #[allow(clippy::needless_pass_by_value)]
    async fn $handler_name(Json(body): Json<Value>) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
      let args: $args_type = serde_json::from_value(body).map_err(json_to_response)?;
      tokio::task::spawn_blocking(move || $api_fn(&args, &NoopProgress))
        .await
        .map_err(join_to_response)?
        .map_err(eyre_to_response)?;
      Ok(Json(serde_json::json!({"status": "ok"})))
=======
struct AppError(Report);

impl IntoResponse for AppError {
  fn into_response(self) -> Response {
    let body = serde_json::json!({
      "code": "internal_error",
      "message": format!("{:?}", self.0),
    });
    (StatusCode::INTERNAL_SERVER_ERROR, Json(body)).into_response()
  }
}

impl<E: Into<Report>> From<E> for AppError {
  fn from(err: E) -> Self {
    Self(err.into())
  }
}

=======
>>>>>>> 759dbc26 (refactor(app-server): extract error, consolidate SSE transport into modules)
async fn handle_version() -> Result<Json<Value>, AppError> {
  let value = serde_json::to_value(version_info())?;
  Ok(Json(value))
}

<<<<<<< HEAD
<<<<<<< HEAD
macro_rules! define_handler {
  ($handler_name:ident, $server_args:ty, $real_args:ty, $api_fn:path, $extract:expr) => {
    #[allow(clippy::needless_pass_by_value)]
    async fn $handler_name(Json(body): Json<Value>) -> Response {
      let args: $server_args = match parse_args(body) {
        Ok(args) => args,
        Err(resp) => return resp,
      };
      let real_args = <$real_args>::from(args);
      sse_response(move |progress| {
        let result = $api_fn(&real_args, progress)?;
        #[allow(clippy::redundant_closure_call)]
        Ok($extract(result))
      })
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)
    }
  };
}

<<<<<<< HEAD
<<<<<<< HEAD
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
=======
define_handler!(handle_ancestral, app_api::TreetimeAncestralArgs, app_api::commands::ancestral);
define_handler!(handle_clock, app_api::TreetimeClockArgs, app_api::commands::clock);
define_handler!(handle_timetree, app_api::TreetimeTimetreeArgs, app_api::commands::timetree);
define_handler!(handle_mugration, app_api::TreetimeMugrationArgs, app_api::commands::mugration);
define_handler!(handle_optimize, app_api::TreetimeOptimizeArgs, app_api::commands::optimize);
define_handler!(handle_prune, app_api::TreetimePruneArgs, app_api::commands::prune);
>>>>>>> 7c23b04d (fix(api,napi,server): clean up Cargo deps, fix error contract, fix code style)
=======
define_handler!(
  handle_ancestral,
  ServerAncestralArgs,
  app_api::TreetimeAncestralArgs,
  app_api::commands::ancestral,
  |result: app_api::AncestralResult| {
    serde_json::json!({
      "model_name": result.model_name.to_string(),
    })
  }
);
=======
async fn handle_ancestral(Json(body): Json<Value>) -> Response {
  handle_command::<ServerAncestralArgs, _, _>(body, app_api::commands::ancestral)
=======
async fn handle_datasets(State(config): State<Arc<ServerConfig>>) -> Result<Json<Value>, AppError> {
  let datasets = discover_datasets(&config.data_dir);
  let value = serde_json::to_value(datasets)?;
  Ok(Json(value))
}

async fn handle_ancestral(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerAncestralArgs, _, _>(body, &config.out_dir, app_api::commands::ancestral)
>>>>>>> 6fc31936 (feat(server): add CLI args for jobs, data dir, and output dir)
}

async fn handle_clock(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerClockArgs, _, _>(body, &config.out_dir, app_api::commands::clock)
}
>>>>>>> 4d100e30 (refactor(app-server): derive Serialize on result types, eliminate handler macro)

async fn handle_timetree(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerTimetreeArgs, _, _>(body, &config.out_dir, app_api::commands::timetree)
}

async fn handle_mugration(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerMugrationArgs, _, _>(body, &config.out_dir, app_api::commands::mugration)
}

async fn handle_optimize(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerOptimizeArgs, _, _>(body, &config.out_dir, app_api::commands::optimize)
}

<<<<<<< HEAD
<<<<<<< HEAD
define_handler!(
  handle_optimize,
  ServerOptimizeArgs,
  app_api::TreetimeOptimizeArgs,
  app_api::commands::optimize,
  |_result: app_api::OptimizeResult| { serde_json::json!({}) }
);

define_handler!(
  handle_prune,
  ServerPruneArgs,
  app_api::TreetimePruneArgs,
  app_api::commands::prune,
  |_result: app_api::PruneResult| { serde_json::json!({}) }
);
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)
=======
async fn handle_prune(Json(body): Json<Value>) -> Response {
  handle_command::<ServerPruneArgs, _, _>(body, app_api::commands::prune)
=======
async fn handle_prune(State(config): State<Arc<ServerConfig>>, Json(body): Json<Value>) -> Response {
  handle_command::<ServerPruneArgs, _, _>(body, &config.out_dir, app_api::commands::prune)
>>>>>>> 6fc31936 (feat(server): add CLI args for jobs, data dir, and output dir)
}
>>>>>>> 4d100e30 (refactor(app-server): derive Serialize on result types, eliminate handler macro)
