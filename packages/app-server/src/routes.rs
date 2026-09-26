use crate::error::AppError;
use crate::events::run_events_sse;
use crate::openapi::{add_component, add_components, add_setting_catalog, schema_ref};
use crate::state::AppState;
use app_commands::bridge::operations::{HttpMethod, HttpOperation, OperationRequest, REQUEST_ARG, http_operations};
use app_commands::config::source::escape_pointer;
use app_commands::job::JobId;
use app_commands::runs::errors::{invalid, parse_request};
use axum::Json;
use axum::Router;
use axum::body::{Body, Bytes};
use axum::extract::{DefaultBodyLimit, Path, Query, RawPathParams, State};
use axum::http::{HeaderMap, HeaderValue, StatusCode, header};
use axum::response::{IntoResponse, Response};
use axum::routing::{MethodFilter, MethodRouter, get, on};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde::Deserialize;
use serde_json::{Map, Value, json};
use std::io;
use std::sync::Arc;
use tokio_stream::StreamExt as _;
use tokio_util::io::{StreamReader, SyncIoBridge};
use treetime_schema::version_info;
use treetime_utils::{make_error, make_report};
use utoipa::openapi::{ContactBuilder, LicenseBuilder};
use utoipa_axum::router::OpenApiRouter;
use utoipa_axum::routes;

const LAST_EVENT_ID: &str = "last-event-id";

const JSON_CONTENT_TYPE: &str = "application/json";

const ERROR_RESPONSES: &str = "default";

const OPERATION_KEY: &str = "x-operation";

const CONTENT_TYPES: &[(&str, &str)] = &[
  ("json", "application/json"),
  ("svg", "image/svg+xml"),
  ("png", "image/png"),
  ("zip", "application/zip"),
  ("nwk", "text/plain; charset=utf-8"),
  ("nexus", "text/plain; charset=utf-8"),
  ("csv", "text/csv; charset=utf-8"),
  ("tsv", "text/tab-separated-values; charset=utf-8"),
  ("fasta", "text/plain; charset=utf-8"),
  ("dot", "text/plain; charset=utf-8"),
  ("jsonl", "application/jsonl"),
];

pub(crate) fn api_routes(state: &Arc<AppState>) -> Router {
  let (router, _api) = api_router().with_state(Arc::clone(state)).split_for_parts();
  let router = http_operations().into_iter().fold(router, |router, operation| {
    let path = operation.path;
    router.route(path, route_operation(operation, Arc::clone(state)))
  });
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

  let mut doc = serde_json::to_value(api)?;
  add_components(&mut doc)?;
  add_setting_catalog(&mut doc)?;
  for operation in http_operations() {
    add_operation(&mut doc, &operation)?;
  }
  for transfer in transfers() {
    describe_transfer(&mut doc, &transfer)?;
  }
  Ok(doc)
}

fn api_router() -> OpenApiRouter<Arc<AppState>> {
  let uploads = OpenApiRouter::new()
    .routes(routes!(handle_upload_input))
    .layer(DefaultBodyLimit::disable());
  OpenApiRouter::new()
    .routes(routes!(handle_operation))
    .routes(routes!(handle_run_events))
    .routes(routes!(handle_run_file))
    .routes(routes!(handle_run_archive))
    .merge(uploads)
}

fn route_operation(operation: HttpOperation, state: Arc<AppState>) -> MethodRouter {
  let filter = match operation.method {
    HttpMethod::Get => MethodFilter::GET,
    HttpMethod::Post => MethodFilter::POST,
    HttpMethod::Patch => MethodFilter::PATCH,
    HttpMethod::Delete => MethodFilter::DELETE,
  };
  let operation = Arc::new(operation);
  on(filter, move |params: RawPathParams, body: Bytes| {
    let operation = Arc::clone(&operation);
    let state = Arc::clone(&state);
    async move { call_operation(&operation, state, &params, &body).await }
  })
}

async fn call_operation(
  operation: &HttpOperation,
  state: Arc<AppState>,
  params: &RawPathParams,
  body: &[u8],
) -> Result<Response, AppError> {
  let mut args: Map<String, Value> = params
    .iter()
    .map(|(name, value)| (name.to_owned(), Value::String(value.to_owned())))
    .collect();
  if operation.body().is_some() {
    let request = if body.is_empty() {
      json!({})
    } else {
      serde_json::from_slice(body).map_err(|err| invalid(err.to_string()))?
    };
    args.insert(REQUEST_ARG.to_owned(), request);
  }
  let request: OperationRequest = parse_request(json!({ "operation": operation.name, "args": args }))?;
  let answer = answer(&state, request).await?;
  Ok(if (operation.response)().is_unit() {
    StatusCode::NO_CONTENT.into_response()
  } else {
    json_response(answer)
  })
}

async fn answer(state: &Arc<AppState>, request: OperationRequest) -> Result<String, Report> {
  let service = Arc::clone(&state.service);
  tokio::task::spawn_blocking(move || request.handle(service.as_ref())).await?
}

#[utoipa::path(
  post,
  path = "/api/operations",
  operation_id = "operationsCall",
  request_body(content = String, content_type = "application/json"),
  responses((status = 200, description = "The JSON result of the operation the request names; `null` for an operation without result"))
)]
async fn handle_operation(State(state): State<Arc<AppState>>, body: Bytes) -> Result<Response, AppError> {
  let request: OperationRequest = serde_json::from_slice(&body).map_err(|err| invalid(err.to_string()))?;
  Ok(json_response(answer(&state, request).await?))
}

fn json_response(json: String) -> Response {
  (
    [(header::CONTENT_TYPE, HeaderValue::from_static(JSON_CONTENT_TYPE))],
    json,
  )
    .into_response()
}

fn add_operation(doc: &mut Value, operation: &HttpOperation) -> Result<(), Report> {
  let parameters = operation
    .path_args()
    .map(|arg| json!({ "name": arg.name, "in": "path", "required": true, "schema": { "type": "string" } }))
    .collect_vec();
  let response = (operation.response)();
  let success = if response.is_unit() {
    json!({ "204": { "description": operation.description } })
  } else {
    add_component(doc, response.name.as_str(), response.schema.clone())?;
    json!({ "200": {
      "description": operation.description,
      "content": { JSON_CONTENT_TYPE: { "schema": schema_ref(&response.name) } },
    }})
  };
  let mut responses = success;
  responses[ERROR_RESPONSES] = json!({
    "description": "The error, with its causes",
    "content": { JSON_CONTENT_TYPE: { "schema": schema_ref("ErrorResponse") } },
  });
  let mut entry = json!({
    "operationId": operation.operation_id,
    "description": operation.description,
    OPERATION_KEY: operation.name,
    "parameters": parameters,
    "responses": responses,
  });
  if let Some(body) = operation.body() {
    let body = (body.schema)();
    add_component(doc, body.name.as_str(), body.schema.clone())?;
    entry["requestBody"] = json!({
      "required": true,
      "content": { JSON_CONTENT_TYPE: { "schema": schema_ref(&body.name) } },
    });
  }
  let method: &'static str = operation.method.into();
  let path_item = doc
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI document must be a JSON object"))?
    .entry("paths")
    .or_insert_with(|| json!({}))
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI paths must be a JSON object"))?
    .entry(operation.path)
    .or_insert_with(|| json!({}))
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI path `{}` must be a JSON object", operation.path))?;
  if path_item.insert(method.to_owned(), entry).is_some() {
    return make_error!("the OpenAPI document has two `{method} {}` operations", operation.path);
  }
  Ok(())
}

struct Transfer {
  path: &'static str,
  method: &'static str,
  request: Option<(&'static str, Value)>,
  response: Option<(&'static str, Value)>,
}

fn transfers() -> Vec<Transfer> {
  let binary = || json!({ "type": "string", "format": "binary" });
  vec![
    Transfer {
      path: "/api/operations",
      method: "post",
      request: Some((JSON_CONTENT_TYPE, schema_ref("OperationRequest"))),
      response: Some((JSON_CONTENT_TYPE, json!({}))),
    },
    Transfer {
      path: "/api/runs/{id}/events",
      method: "get",
      request: None,
      response: Some(("text/event-stream", schema_ref("RunEvent"))),
    },
    Transfer {
      path: "/api/runs/{id}/inputs/{name}",
      method: "put",
      request: Some(("application/octet-stream", binary())),
      response: Some((JSON_CONTENT_TYPE, schema_ref("UploadedInput"))),
    },
    Transfer {
      path: "/api/runs/{id}/file",
      method: "get",
      request: None,
      response: Some(("application/octet-stream", binary())),
    },
    Transfer {
      path: "/api/runs/{id}/archive",
      method: "get",
      request: None,
      response: Some(("application/zip", binary())),
    },
  ]
}

fn describe_transfer(doc: &mut Value, transfer: &Transfer) -> Result<(), Report> {
  let pointer = format!("/paths/{}/{}", escape_pointer(transfer.path), transfer.method);
  let entry = doc
    .pointer_mut(&pointer)
    .and_then(Value::as_object_mut)
    .ok_or_else(|| {
      make_report!(
        "the OpenAPI document has no operation `{} {}`",
        transfer.method,
        transfer.path
      )
    })?;
  if let Some(responses) = entry.get_mut("responses").and_then(Value::as_object_mut) {
    for (status, response) in responses.iter_mut() {
      if status.starts_with('4') || status.starts_with('5') {
        if let Some(response) = response.as_object_mut() {
          response.insert(
            "content".to_owned(),
            json!({ JSON_CONTENT_TYPE: { "schema": schema_ref("ErrorResponse") } }),
          );
        }
      }
    }
  }
  if let Some((content_type, schema)) = &transfer.request {
    let mut content = Map::new();
    content.insert((*content_type).to_owned(), json!({ "schema": schema }));
    entry.insert(
      "requestBody".to_owned(),
      json!({ "required": true, "content": Value::Object(content) }),
    );
  }
  if let Some((content_type, schema)) = &transfer.response {
    let response = entry
      .get_mut("responses")
      .and_then(|responses| responses.get_mut("200"))
      .and_then(Value::as_object_mut)
      .ok_or_else(|| {
        make_report!(
          "operation `{} {}` has no success response",
          transfer.method,
          transfer.path
        )
      })?;
    let mut content = Map::new();
    content.insert((*content_type).to_owned(), json!({ "schema": schema }));
    response.insert("content".to_owned(), Value::Object(content));
  }
  Ok(())
}

#[derive(Deserialize)]
struct EventsQuery {
  from: Option<usize>,
}

#[utoipa::path(
  get,
  path = "/api/runs/{id}/events",
  operation_id = "runsEvents",
  params(
    ("id" = String, Path, description = "Id of the run"),
    ("from" = Option<usize>, Query, description = "Sequence number of the first event to send; the `Last-Event-ID` header of a reconnect resumes after that event"),
  ),
  responses((status = 200, description = "Stream of the run's events from `from`: `started`, `progress`, `log` and `iteration`, then one `terminal`"))
)]
async fn handle_run_events(
  State(state): State<Arc<AppState>>,
  Path(id): Path<String>,
  Query(query): Query<EventsQuery>,
  headers: HeaderMap,
) -> Result<Response, AppError> {
  let resume = headers
    .get(LAST_EVENT_ID)
    .and_then(|value| value.to_str().ok())
    .and_then(|value| value.parse::<usize>().ok())
    .map(|seq| seq + 1);
  let from = resume.or(query.from).unwrap_or(0);
  run_events_sse(&state, &run_id(&id)?, from)
}

#[utoipa::path(
  put,
  path = "/api/runs/{id}/inputs/{name}",
  operation_id = "runsUploadInput",
  params(
    ("id" = String, Path, description = "Id of a run that has not started"),
    ("name" = String, Path, description = "File name inside the run's `inputs/` folder"),
  ),
  responses(
    (status = 200, description = "The stored file and the path to use for it in the run's configuration"),
    (status = 413, description = "The inputs of the run exceed the upload limit of the server"),
  )
)]
async fn handle_upload_input(
  State(state): State<Arc<AppState>>,
  Path((id, name)): Path<(String, String)>,
  body: Body,
) -> Result<Json<Value>, AppError> {
  let id = run_id(&id)?;
  let stream = body.into_data_stream().map(|chunk| chunk.map_err(io::Error::other));
  let mut reader = SyncIoBridge::new(StreamReader::new(stream));
  let runs = Arc::clone(&state.runs);
  let limit = state.config.max_upload_size;
  let uploaded = tokio::task::spawn_blocking(move || runs.upload_input(&id, &name, &mut reader, limit)).await??;
  Ok(Json(serde_json::to_value(uploaded)?))
}

#[derive(Deserialize)]
struct FileQuery {
  path: String,
}

#[utoipa::path(
  get,
  path = "/api/runs/{id}/file",
  operation_id = "runsFile",
  params(
    ("id" = String, Path, description = "Id of the run"),
    ("path" = String, Query, description = "Path of the file relative to the run's `out/` folder"),
  ),
  responses((status = 200, description = "Contents of the file"))
)]
async fn handle_run_file(
  State(state): State<Arc<AppState>>,
  Path(id): Path<String>,
  Query(query): Query<FileQuery>,
) -> Result<Response, AppError> {
  let path = state.runs.file_path(&run_id(&id)?, &query.path)?;
  let bytes = tokio::fs::read(&path)
    .await
    .wrap_err_with(|| format!("When reading '{}'", path.display()))?;
  let extension = path.extension().and_then(|extension| extension.to_str()).unwrap_or("");
  Ok(([(header::CONTENT_TYPE, content_type(extension))], bytes).into_response())
}

#[utoipa::path(
  get,
  path = "/api/runs/{id}/archive",
  operation_id = "runsArchive",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 200, description = "Zip archive of the run's `out/` folder"))
)]
async fn handle_run_archive(State(state): State<Arc<AppState>>, Path(id): Path<String>) -> Result<Response, AppError> {
  let id = run_id(&id)?;
  let runs = Arc::clone(&state.runs);
  let archive_id = id.clone();
  let bytes = tokio::task::spawn_blocking(move || runs.zip(&archive_id)).await??;
  let disposition = HeaderValue::from_str(&format!("attachment; filename=\"treetime-{}.zip\"", id.as_str()))?;
  Ok(
    (
      [
        (header::CONTENT_TYPE, HeaderValue::from_static("application/zip")),
        (header::CONTENT_DISPOSITION, disposition),
      ],
      bytes,
    )
      .into_response(),
  )
}

fn run_id(id: &str) -> Result<JobId, Report> {
  JobId::parse(id).map_err(|err| invalid(err.to_string()))
}

fn content_type(extension: &str) -> HeaderValue {
  let content_type = CONTENT_TYPES
    .iter()
    .find(|(known, _)| *known == extension)
    .map_or("application/octet-stream", |(_, content_type)| *content_type);
  HeaderValue::from_static(content_type)
}
