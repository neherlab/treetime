<<<<<<< HEAD
<<<<<<< HEAD
use app_api::progress::NoopProgress;
<<<<<<< HEAD
=======
use crate::args::{
  ServerAncestralArgs, ServerClockArgs, ServerMugrationArgs, ServerOptimizeArgs, ServerPruneArgs, ServerTimetreeArgs,
};
use app_api::datasets::discover_datasets;
use crate::error::AppError;
use crate::sse::handle_command;
use crate::state::ServerConfig;
use app_api::datasets::discover_datasets;
=======
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
use app_api::version::version_info;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> f8a8231c (feat(app-server): wire real computation with channel-based SSE progress)
=======
use app_api::version::version_info;
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
use axum::http::StatusCode;
use axum::response::sse::{Event, Sse};
use axum::response::{IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
<<<<<<< HEAD
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
=======
use eyre::Report;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use smart_default::SmartDefault;
use std::convert::Infallible;
use std::time::Duration;
use tokio_stream::StreamExt as _;
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
=======
use std::path::Path;
>>>>>>> 3b85da3a (feat(app): add dataset discovery endpoint)

pub fn api_routes() -> Router {
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
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

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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

<<<<<<< HEAD
#[derive(Serialize)]
struct ProgressEvent {
  stage: String,
  fraction: f64,
  message: String,
}

fn sse_response(phases: &'static [(&str, f64)], result: Value) -> Response {
  let stream = tokio_stream::iter(phases.iter().copied())
    .then(|(stage, fraction)| async move {
      tokio::time::sleep(Duration::from_millis(400)).await;
      let payload = ProgressEvent {
        stage: stage.to_owned(),
        fraction,
        message: stage.to_owned(),
      };
      Ok::<_, Infallible>(
        Event::default()
          .event("progress")
          .json_data(payload)
          .expect("serialization"),
      )
    })
    .chain(tokio_stream::once(Ok::<_, Infallible>(
      Event::default()
        .event("result")
        .json_data(result)
        .expect("serialization"),
    )));

  Sse::new(stream).into_response()
}

// TODO: replace with real args types once serde(default) is added to TreetimeArgs structs

#[derive(SmartDefault, Deserialize)]
#[serde(default)]
struct AncestralArgs {
  tree: String,
  outdir: String,
  input_fastas: Vec<String>,
  aln: Option<String>,
  #[default = "marginal"]
  method_anc: String,
  #[default = "infer"]
  model_name: String,
  #[default = "nuc"]
  alphabet: String,
  dense: Option<bool>,
  aa: Option<bool>,
}

const ANCESTRAL_PHASES: &[(&str, f64)] = &[
  ("Reading alignment", 0.1),
  ("Parsing tree", 0.2),
  ("Inferring GTR model", 0.4),
  ("Running Fitch parsimony", 0.6),
  ("Marginal reconstruction", 0.8),
  ("Writing output", 1.0),
];

fn run_ancestral(_args: &AncestralArgs) -> eyre::Result<Value> {
  // TODO: convert to TreetimeAncestralArgs, call run_ancestral_reconstruction
  Ok(serde_json::json!({
    "model_name": "GTR (inferred)",
    "gtr": {
      "pi": [0.334, 0.198, 0.223, 0.245],
      "W": [[0.0, 0.94, 2.41, 0.52], [0.94, 0.0, 0.48, 2.68], [2.41, 0.48, 0.0, 0.87], [0.52, 2.68, 0.87, 0.0]],
      "mu": 1.0
    },
    "n_sequences": 20,
    "n_sites": 1701
  }))
}

#[derive(SmartDefault, Deserialize)]
#[serde(default)]
struct ClockArgs {
  tree: String,
  outdir: String,
  dates: String,
  aln: Vec<String>,
  #[default = "least-squares"]
  reroot: String,
  #[default = 3.0]
  clock_filter: f64,
  keep_root: Option<bool>,
  covariation: Option<bool>,
}

const CLOCK_PHASES: &[(&str, f64)] = &[
  ("Reading metadata", 0.1),
  ("Parsing tree", 0.2),
  ("Computing root-to-tip distances", 0.4),
  ("Linear regression", 0.6),
  ("Rerooting (least squares)", 0.8),
  ("Writing output", 1.0),
];

fn run_clock(_args: &ClockArgs) -> eyre::Result<Value> {
  // TODO: convert to TreetimeClockArgs, call run_clock
  Ok(serde_json::json!({
    "clock_model": {
      "slope": 3.3e-3,
      "intercept": -6.57,
      "r_squared": 0.97,
      "chi_squared": 8.7,
      "root_date": 1998.7
    },
    "regression_results": [
      { "name": "A/New_York/182/2000", "date": 2000.134, "divergence": 0.0021, "predicted": 0.0044, "outlier": false },
      { "name": "A/Canterbury/58/2000", "date": 2000.682, "divergence": 0.0072, "predicted": 0.0062, "outlier": false },
      { "name": "A/Denmark/107/2003", "date": 2003.003, "divergence": 0.0253, "predicted": 0.0139, "outlier": true },
      { "name": "A/Hawaii/02/2013", "date": 2013.405, "divergence": 0.0441, "predicted": 0.0482, "outlier": false }
    ],
    "n_tips": 20,
    "n_outliers": 2
  }))
}

#[derive(SmartDefault, Deserialize)]
#[serde(default)]
struct TimetreeArgs {
  tree: String,
  outdir: String,
  dates: String,
  input_fastas: Vec<String>,
  aln: Option<String>,
  clock_rate: Option<f64>,
  #[default = 3]
  max_iter: usize,
}

const TIMETREE_PHASES: &[(&str, f64)] = &[
  ("Reading alignment and metadata", 0.05),
  ("Parsing tree", 0.1),
  ("Ancestral reconstruction", 0.2),
  ("Clock regression", 0.3),
  ("Rerooting", 0.4),
  ("Estimating branch lengths (iter 1/3)", 0.55),
  ("Estimating branch lengths (iter 2/3)", 0.7),
  ("Estimating branch lengths (iter 3/3)", 0.85),
  ("Confidence intervals", 0.95),
  ("Writing output", 1.0),
];

fn run_timetree(_args: &TimetreeArgs) -> eyre::Result<Value> {
  // TODO: convert to TreetimeTimetreeArgs, call run_timetree_estimation
  Ok(serde_json::json!({
    "clock_model": {
      "slope": 3.3e-3,
      "intercept": -6.57,
      "r_squared": 0.97
    },
    "n_iterations": 3,
    "log_likelihood": -12345.6,
    "tmrca": 1998.7
  }))
}

#[derive(SmartDefault, Deserialize)]
#[serde(default)]
struct MugrationArgs {
  tree: String,
  outdir: String,
  states: String,
  #[default = "country"]
  attribute: String,
  #[default = "strain"]
  name_column: String,
  #[default = 1.0]
  pc: f64,
  missing_data: Option<String>,
}

const MUGRATION_PHASES: &[(&str, f64)] = &[
  ("Reading metadata", 0.1),
  ("Parsing tree", 0.2),
  ("Inferring transition model", 0.4),
  ("Marginal reconstruction", 0.7),
  ("Computing confidences", 0.9),
  ("Writing output", 1.0),
];

fn run_mugration(_args: &MugrationArgs) -> eyre::Result<Value> {
  // TODO: convert to TreetimeMugrationArgs, call run_mugration
  Ok(serde_json::json!({
    "log_lh": -42.3,
    "n_states": 8,
    "attribute": "country",
    "alphabet": ["USA", "Hong Kong", "Denmark", "New Zealand", "Nicaragua", "Viet Nam", "Peru", "Iran"]
  }))
}

#[derive(SmartDefault, Deserialize)]
#[serde(default)]
struct OptimizeArgs {
  tree: String,
  outdir: String,
  input_fastas: Vec<String>,
  aln: Option<String>,
  #[default = "infer"]
  model_name: String,
  #[default = 10]
  max_iter: usize,
}

const OPTIMIZE_PHASES: &[(&str, f64)] = &[
  ("Reading alignment", 0.1),
  ("Parsing tree", 0.2),
  ("Initial likelihood", 0.3),
  ("Optimizing branch lengths (iter 1)", 0.5),
  ("Optimizing branch lengths (iter 2)", 0.7),
  ("Optimizing branch lengths (iter 3)", 0.9),
  ("Writing output", 1.0),
];

fn run_optimize(_args: &OptimizeArgs) -> eyre::Result<Value> {
  // TODO: convert to TreetimeOptimizeArgs, call run_optimize
  Ok(serde_json::json!({
    "n_iterations": 3,
    "initial_log_lh": -13500.0,
    "final_log_lh": -12345.6,
    "improvement": 1154.4
  }))
}

#[derive(SmartDefault, Deserialize)]
#[serde(default)]
struct PruneArgs {
  tree: String,
  outdir: String,
  input_fastas: Vec<String>,
  aln: Option<String>,
  prune_short: Option<f64>,
  prune_empty: bool,
}

const PRUNE_PHASES: &[(&str, f64)] = &[
  ("Reading alignment", 0.15),
  ("Parsing tree", 0.3),
  ("Identifying short branches", 0.5),
  ("Pruning nodes", 0.75),
  ("Writing output", 1.0),
];

fn run_prune(_args: &PruneArgs) -> eyre::Result<Value> {
  // TODO: convert to TreetimePruneArgs, call run_prune
  Ok(serde_json::json!({
    "nodes_removed": 4,
    "branches_removed": 4,
    "nodes_before": 37,
    "nodes_after": 33,
    "leaves_before": 19,
    "leaves_after": 17
  }))
}

>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
macro_rules! define_handler {
  ($handler_name:ident, $args_type:ty, $phases:expr, $run_fn:path) => {
    #[allow(clippy::needless_pass_by_value)]
<<<<<<< HEAD
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
=======
    async fn $handler_name(Json(body): Json<Value>) -> Response {
      let args: $args_type = match serde_json::from_value(body) {
        Ok(args) => args,
        Err(err) => return AppError::from(err).into_response(),
      };
      match $run_fn(&args) {
        Ok(result) => sse_response($phases, result),
        Err(err) => AppError(err).into_response(),
      }
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
    }
  };
}

<<<<<<< HEAD
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
=======
async fn handle_datasets() -> Result<Json<Value>, AppError> {
  let data_dir = std::env::var("DATA_DIR").unwrap_or_else(|_| "data".to_owned());
  let datasets = discover_datasets(Path::new(&data_dir));
  let value = serde_json::to_value(datasets)?;
  Ok(Json(value))
}

>>>>>>> 3b85da3a (feat(app): add dataset discovery endpoint)
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
=======
define_handler!(handle_ancestral, AncestralArgs, ANCESTRAL_PHASES, run_ancestral);
define_handler!(handle_clock, ClockArgs, CLOCK_PHASES, run_clock);
define_handler!(handle_timetree, TimetreeArgs, TIMETREE_PHASES, run_timetree);
define_handler!(handle_mugration, MugrationArgs, MUGRATION_PHASES, run_mugration);
define_handler!(handle_optimize, OptimizeArgs, OPTIMIZE_PHASES, run_optimize);
define_handler!(handle_prune, PruneArgs, PRUNE_PHASES, run_prune);
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
