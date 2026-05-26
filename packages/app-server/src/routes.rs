<<<<<<< HEAD
use app_api::progress::NoopProgress;
=======
use crate::args::{
  ServerAncestralArgs, ServerClockArgs, ServerMugrationArgs, ServerOptimizeArgs, ServerPruneArgs, ServerTimetreeArgs,
};
use crate::progress::ChannelProgress;
use app_api::progress::ProgressSink;
use app_api::version::version_info;
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)
use axum::http::StatusCode;
use axum::routing::post;
use axum::{Json, Router};
<<<<<<< HEAD
use serde_json::Value;

fn error_response(status: StatusCode, message: String) -> (StatusCode, Json<Value>) {
  (status, Json(serde_json::json!({"status": "error", "error": message})))
}

fn eyre_to_response(err: eyre::Report) -> (StatusCode, Json<Value>) {
  error_response(StatusCode::INTERNAL_SERVER_ERROR, format!("{err:?}"))
}

fn json_to_response(err: serde_json::Error) -> (StatusCode, Json<Value>) {
  error_response(StatusCode::BAD_REQUEST, format!("{err}"))
}

fn join_to_response(err: tokio::task::JoinError) -> (StatusCode, Json<Value>) {
  error_response(StatusCode::INTERNAL_SERVER_ERROR, format!("{err}"))
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

pub fn api_routes() -> Router {
  Router::new()
    .route("/ancestral", post(handle_ancestral))
    .route("/clock", post(handle_clock))
    .route("/timetree", post(handle_timetree))
    .route("/mugration", post(handle_mugration))
    .route("/optimize", post(handle_optimize))
    .route("/prune", post(handle_prune))
}

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

async fn handle_version() -> Result<Json<Value>, AppError> {
  let value = serde_json::to_value(version_info())?;
  Ok(Json(value))
}

#[allow(tail_expr_drop_order)]
fn sse_response<F>(run_fn: F) -> Response
where
  F: FnOnce(&dyn ProgressSink) -> Result<Value, Report> + Send + 'static,
{
  let (tx, rx) = mpsc::unbounded_channel::<ProgressEvent>();

  let computation = tokio::task::spawn_blocking(move || {
    let progress = ChannelProgress::new(tx);
    run_fn(&progress)
  });

  let stream = async_stream::stream! {
    let mut rx_stream = UnboundedReceiverStream::new(rx);
    while let Some(event) = rx_stream.next().await {
      yield Ok::<_, Infallible>(
        Event::default()
          .event("progress")
          .json_data(event)
          .expect("ProgressEvent serialization"),
      );
    }

    let result_value = match computation.await {
      Ok(Ok(value)) => value,
      Ok(Err(err)) => serde_json::json!({ "error": format!("{err:?}") }),
      Err(err) => serde_json::json!({ "error": format!("computation panicked: {err}") }),
    };
    yield Ok::<_, Infallible>(
      Event::default()
        .event("result")
        .json_data(result_value)
        .expect("result serialization"),
    );
  };

  Sse::new(stream).into_response()
}

#[allow(clippy::result_large_err)]
fn parse_args<A: DeserializeOwned>(body: Value) -> Result<A, Response> {
  serde_json::from_value(body).map_err(|err| AppError::from(err).into_response())
}

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

define_handler!(
  handle_clock,
  ServerClockArgs,
  app_api::TreetimeClockArgs,
  app_api::commands::clock,
  |result: app_api::ClockResult| {
    serde_json::json!({
      "clock_model": {
        "slope": result.clock_model.clock_rate(),
        "intercept": result.clock_model.intercept(),
      },
    })
  }
);

define_handler!(
  handle_timetree,
  ServerTimetreeArgs,
  app_api::TreetimeTimetreeArgs,
  app_api::commands::timetree,
  |result: app_api::TimetreeResult| {
    serde_json::json!({
      "clock_model": {
        "slope": result.clock_model.clock_rate(),
        "intercept": result.clock_model.intercept(),
      },
    })
  }
);

define_handler!(
  handle_mugration,
  ServerMugrationArgs,
  app_api::TreetimeMugrationArgs,
  app_api::commands::mugration,
  |result: app_api::MugrationResult| {
    serde_json::json!({
      "log_lh": result.log_lh,
    })
  }
);

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
