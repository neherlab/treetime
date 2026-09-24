use crate::contract::{CancelJobResponse, DatasetInfo, ErrorResponse};
use crate::error::AppError;
use crate::openapi::{add_components, config_component, schema_ref};
use crate::sse::run_command_sse;
use crate::state::AppState;
use app_commands::command::{AppCommand, CheckConfigRequest, check_config};
use app_commands::job::JobId;
use app_datasets::discover_datasets;
use axum::extract::{Path, State};
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::routing::get;
use axum::{Json, Router};
use eyre::Report;
use serde_json::{Map, Value, json};
use std::sync::Arc;
use strum::IntoEnumIterator;
use treetime_schema::version_info;
use treetime_utils::make_report;
use utoipa::OpenApi;
use utoipa::openapi::{ContactBuilder, LicenseBuilder};
use utoipa_axum::router::OpenApiRouter;
use utoipa_axum::routes;

pub(crate) fn api_routes(state: Arc<AppState>) -> Router {
  let (router, _api) = api_router().with_state(state).split_for_parts();
  router
    .route("/api/health", get(handle_health))
    .route("/api/openapi.json", get(handle_openapi))
}

async fn handle_health() -> Json<Value> {
  Json(json!({
    "status": "ok",
    "version": version_info().version,
  }))
}

async fn handle_openapi() -> Result<Json<Value>, AppError> {
  Ok(Json(api_doc()?))
}

pub fn api_doc() -> Result<Value, Report> {
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

  let mut doc = serde_json::to_value(api)?;
  add_components(&mut doc)?;
  for operation in operations() {
    describe_operation(&mut doc, &operation)?;
  }
  Ok(doc)
}

fn api_router() -> OpenApiRouter<Arc<AppState>> {
  OpenApiRouter::new()
    .routes(routes!(handle_version))
    .routes(routes!(handle_datasets))
    .routes(routes!(handle_check_config))
    .routes(routes!(handle_cancel_job))
    .routes(routes!(handle_timetree))
    .routes(routes!(handle_optimize))
    .routes(routes!(handle_prune))
    .routes(routes!(handle_ancestral))
    .routes(routes!(handle_clock))
    .routes(routes!(handle_mugration))
}

struct Operation {
  path: String,
  method: &'static str,
  bridge_type: &'static str,
  request: Option<Value>,
  response: Option<(&'static str, Value)>,
}

fn operations() -> Vec<Operation> {
  let fixed = [
    Operation {
      path: "/api/version".to_owned(),
      method: "get",
      bridge_type: "query",
      request: None,
      response: Some(("application/json", schema_ref("VersionInfo"))),
    },
    Operation {
      path: "/api/datasets".to_owned(),
      method: "get",
      bridge_type: "query",
      request: None,
      response: None,
    },
    Operation {
      path: "/api/check-config".to_owned(),
      method: "post",
      bridge_type: "request",
      request: Some(schema_ref("CheckConfigRequest")),
      response: Some(("application/json", schema_ref("CheckConfigResponse"))),
    },
    Operation {
      path: "/api/jobs/{job_id}/cancel".to_owned(),
      method: "post",
      bridge_type: "request",
      request: None,
      response: None,
    },
  ];
  let commands = AppCommand::iter().map(|command| Operation {
    path: format!("/api/{command}"),
    method: "post",
    bridge_type: "command",
    request: Some(schema_ref(&config_component(command))),
    response: Some(("text/event-stream", schema_ref("JobEvent"))),
  });
  fixed.into_iter().chain(commands).collect()
}

fn describe_operation(doc: &mut Value, operation: &Operation) -> Result<(), Report> {
  let pointer = format!(
    "/paths/{}/{}",
    operation.path.replace('~', "~0").replace('/', "~1"),
    operation.method
  );
  let entry = doc
    .pointer_mut(&pointer)
    .and_then(Value::as_object_mut)
    .ok_or_else(|| {
      make_report!(
        "the OpenAPI document has no operation `{} {}`",
        operation.method,
        operation.path
      )
    })?;
  entry.insert("x-bridge-type".to_owned(), json!(operation.bridge_type));
  if let Some(schema) = &operation.request {
    entry.insert(
      "requestBody".to_owned(),
      json!({ "required": true, "content": { "application/json": { "schema": schema } } }),
    );
  }
  if let Some((content_type, schema)) = &operation.response {
    let response = entry
      .get_mut("responses")
      .and_then(|responses| responses.get_mut("200"))
      .and_then(Value::as_object_mut)
      .ok_or_else(|| make_report!("operation `{}` has no success response", operation.path))?;
    let mut content = Map::new();
    content.insert((*content_type).to_owned(), json!({ "schema": schema }));
    response.insert("content".to_owned(), Value::Object(content));
  }
  Ok(())
}

#[utoipa::path(
  get,
  path = "/api/version",
  operation_id = "version",
  responses((status = 200, description = "Version information"))
)]
async fn handle_version() -> Result<Json<Value>, AppError> {
  Ok(Json(serde_json::to_value(version_info())?))
}

#[utoipa::path(
  get,
  path = "/api/datasets",
  operation_id = "datasets",
  responses((status = 200, description = "Example datasets in the data directory", body = Vec<DatasetInfo>))
)]
async fn handle_datasets(State(state): State<Arc<AppState>>) -> Result<Json<Value>, AppError> {
  let datasets = discover_datasets(&state.config.data_dir)?;
  Ok(Json(serde_json::to_value(datasets)?))
}

#[utoipa::path(
  post,
  path = "/api/check-config",
  operation_id = "configCheck",
  responses((status = 200, description = "The configuration with every default filled in, or the problems found in it"))
)]
async fn handle_check_config(Json(body): Json<Value>) -> Response {
  match serde_json::from_value::<CheckConfigRequest>(body) {
    Ok(request) => Json(check_config(&request)).into_response(),
    Err(err) => error_response(StatusCode::BAD_REQUEST, "invalid_request", &Report::from(err)),
  }
}

#[utoipa::path(
  post,
  path = "/api/jobs/{job_id}/cancel",
  operation_id = "jobCancel",
  params(("job_id" = String, Path, description = "Id of the job, from its `started` event")),
  responses(
    (status = 200, description = "Cancellation requested; the job ends with a `cancelled` terminal event", body = CancelJobResponse),
    (status = 404, description = "No job with this id is running", body = ErrorResponse),
  )
)]
async fn handle_cancel_job(State(state): State<Arc<AppState>>, Path(job_id): Path<String>) -> Response {
  let job_id = match JobId::parse(&job_id) {
    Ok(job_id) => job_id,
    Err(err) => return error_response(StatusCode::BAD_REQUEST, "invalid_job_id", &err),
  };
  if state.jobs.cancel(&job_id) {
    Json(CancelJobResponse { cancelled: true }).into_response()
  } else {
    error_response(
      StatusCode::NOT_FOUND,
      "job_not_found",
      &make_report!("no job with id `{}` is running", job_id.as_str()),
    )
  }
}

fn error_response(status: StatusCode, code: &str, err: &Report) -> Response {
  (status, Json(json!({ "code": code, "message": err.to_string() }))).into_response()
}

macro_rules! command_route {
  ($handler:ident, $command:ident, $path:literal, $operation_id:literal) => {
    #[utoipa::path(
      post,
      path = $path,
      operation_id = $operation_id,
      responses((status = 200, description = "Stream of job events: `started`, then `progress` and `log`, then one `terminal`"))
    )]
    async fn $handler(State(state): State<Arc<AppState>>, Json(config): Json<Value>) -> Response {
      run_command_sse(&state, AppCommand::$command, config)
    }
  };
}

command_route!(handle_timetree, Timetree, "/api/timetree", "timetree");
command_route!(handle_optimize, Optimize, "/api/optimize", "optimize");
command_route!(handle_prune, Prune, "/api/prune", "prune");
command_route!(handle_ancestral, Ancestral, "/api/ancestral", "ancestral");
command_route!(handle_clock, Clock, "/api/clock", "clock");
command_route!(handle_mugration, Mugration, "/api/mugration", "mugration");

#[derive(OpenApi)]
#[openapi(components(schemas(ErrorResponse)))]
struct SharedSchemas;
