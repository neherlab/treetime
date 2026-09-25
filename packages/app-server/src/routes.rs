use crate::contract::{CancelRunResponse, ErrorResponse};
use crate::error::AppError;
use crate::events::run_events_sse;
use crate::openapi::{add_components, add_setting_catalog, schema_ref};
use crate::state::AppState;
use app_commands::check_config::{CheckConfigRequest, check_config};
use app_commands::check_inputs::{CheckInputsRequest, check_inputs};
use app_commands::command::AppCommand;
use app_commands::job::JobId;
use app_commands::run_config::{RunConfigRequest, run_config};
use app_commands::runs::errors::invalid;
use app_commands::runs::record::{CreateRunRequest, RunRecord, StartRunRequest, UpdateRunRequest};
use app_datasets::discover_datasets;
use axum::body::Body;
use axum::extract::{DefaultBodyLimit, Path, Query, State};
use axum::http::{HeaderMap, HeaderValue, StatusCode, header};
use axum::response::{IntoResponse, Response};
use axum::routing::get;
use axum::{Json, Router};
use eyre::{Report, WrapErr};
use log::{error, info};
use serde::Deserialize;
use serde_json::{Map, Value, json};
use std::io;
use std::path::Path as FsPath;
use std::sync::Arc;
use strum::VariantNames;
use tokio_stream::StreamExt as _;
use tokio_util::io::{StreamReader, SyncIoBridge};
use treetime_schema::version_info;
use treetime_utils::make_report;
use utoipa::OpenApi;
use utoipa::openapi::{ContactBuilder, LicenseBuilder};
use utoipa_axum::router::OpenApiRouter;
use utoipa_axum::routes;

const LAST_EVENT_ID: &str = "last-event-id";

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
  add_setting_catalog(&mut doc)?;
  for operation in operations() {
    describe_operation(&mut doc, &operation)?;
  }
  Ok(doc)
}

fn api_router() -> OpenApiRouter<Arc<AppState>> {
  let uploads = OpenApiRouter::new()
    .routes(routes!(handle_upload_input))
    .layer(DefaultBodyLimit::disable());
  OpenApiRouter::new()
    .routes(routes!(handle_version))
    .routes(routes!(handle_datasets))
    .routes(routes!(handle_check_config))
    .routes(routes!(handle_run_config))
    .routes(routes!(handle_check_inputs))
    .routes(routes!(handle_list_runs, handle_create_run))
    .routes(routes!(handle_get_run, handle_update_run, handle_delete_run))
    .routes(routes!(handle_start_run))
    .routes(routes!(handle_cancel_run))
    .routes(routes!(handle_restore_run))
    .routes(routes!(handle_purge_run))
    .routes(routes!(handle_run_events))
    .routes(routes!(handle_run_files))
    .routes(routes!(handle_run_file))
    .routes(routes!(handle_run_archive))
    .merge(uploads)
}

struct Operation {
  path: &'static str,
  method: &'static str,
  bridge_type: &'static str,
  request: Option<(&'static str, Value)>,
  response: Option<(&'static str, Value)>,
}

fn operations() -> Vec<Operation> {
  let json_body = |component: &str| Some(("application/json", schema_ref(component)));
  vec![
    Operation {
      path: "/api/version",
      method: "get",
      bridge_type: "query",
      request: None,
      response: json_body("VersionInfo"),
    },
    Operation {
      path: "/api/datasets",
      method: "get",
      bridge_type: "query",
      request: None,
      response: json_body("DatasetCatalog"),
    },
    Operation {
      path: "/api/check-config",
      method: "post",
      bridge_type: "request",
      request: json_body("CheckConfigRequest"),
      response: json_body("CheckConfigResponse"),
    },
    Operation {
      path: "/api/run-config",
      method: "post",
      bridge_type: "request",
      request: json_body("RunConfigRequest"),
      response: json_body("RunConfigResponse"),
    },
    Operation {
      path: "/api/check-inputs",
      method: "post",
      bridge_type: "request",
      request: json_body("CheckInputsRequest"),
      response: json_body("InputFacts"),
    },
    Operation {
      path: "/api/runs",
      method: "get",
      bridge_type: "query",
      request: None,
      response: json_body("RunList"),
    },
    Operation {
      path: "/api/runs",
      method: "post",
      bridge_type: "request",
      request: json_body("CreateRunRequest"),
      response: json_body("RunRecord"),
    },
    Operation {
      path: "/api/runs/{id}",
      method: "get",
      bridge_type: "query",
      request: None,
      response: json_body("RunRecord"),
    },
    Operation {
      path: "/api/runs/{id}",
      method: "patch",
      bridge_type: "request",
      request: json_body("UpdateRunRequest"),
      response: json_body("RunSummary"),
    },
    Operation {
      path: "/api/runs/{id}",
      method: "delete",
      bridge_type: "request",
      request: None,
      response: None,
    },
    Operation {
      path: "/api/runs/{id}/start",
      method: "post",
      bridge_type: "request",
      request: json_body("StartRunRequest"),
      response: json_body("RunRecord"),
    },
    Operation {
      path: "/api/runs/{id}/cancel",
      method: "post",
      bridge_type: "request",
      request: None,
      response: json_body("CancelRunResponse"),
    },
    Operation {
      path: "/api/runs/{id}/restore",
      method: "post",
      bridge_type: "request",
      request: None,
      response: json_body("RunSummary"),
    },
    Operation {
      path: "/api/runs/{id}/purge",
      method: "post",
      bridge_type: "request",
      request: None,
      response: None,
    },
    Operation {
      path: "/api/runs/{id}/events",
      method: "get",
      bridge_type: "stream",
      request: None,
      response: Some(("text/event-stream", schema_ref("RunEvent"))),
    },
    Operation {
      path: "/api/runs/{id}/inputs/{name}",
      method: "put",
      bridge_type: "upload",
      request: Some((
        "application/octet-stream",
        json!({ "type": "string", "format": "binary" }),
      )),
      response: json_body("UploadedInput"),
    },
    Operation {
      path: "/api/runs/{id}/files",
      method: "get",
      bridge_type: "query",
      request: None,
      response: Some((
        "application/json",
        json!({ "type": "array", "items": schema_ref("RunFile") }),
      )),
    },
    Operation {
      path: "/api/runs/{id}/file",
      method: "get",
      bridge_type: "file",
      request: None,
      response: Some((
        "application/octet-stream",
        json!({ "type": "string", "format": "binary" }),
      )),
    },
    Operation {
      path: "/api/runs/{id}/archive",
      method: "get",
      bridge_type: "file",
      request: None,
      response: Some(("application/zip", json!({ "type": "string", "format": "binary" }))),
    },
  ]
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
  if let Some((content_type, schema)) = &operation.request {
    let mut content = Map::new();
    content.insert((*content_type).to_owned(), json!({ "schema": schema }));
    entry.insert(
      "requestBody".to_owned(),
      json!({ "required": true, "content": Value::Object(content) }),
    );
  }
  if let Some((content_type, schema)) = &operation.response {
    let response = entry
      .get_mut("responses")
      .and_then(|responses| responses.get_mut("200"))
      .and_then(Value::as_object_mut)
      .ok_or_else(|| {
        make_report!(
          "operation `{} {}` has no success response",
          operation.method,
          operation.path
        )
      })?;
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
  responses((status = 200, description = "Example datasets and example configurations in the data directory"))
)]
async fn handle_datasets(State(state): State<Arc<AppState>>) -> Result<Json<Value>, AppError> {
  let catalog = discover_datasets(&state.config.data_dir, AppCommand::VARIANTS)?;
  Ok(Json(serde_json::to_value(catalog)?))
}

#[utoipa::path(
  post,
  path = "/api/check-config",
  operation_id = "configCheck",
  responses((status = 200, description = "The configuration with every default filled in, or the problems found in it"))
)]
async fn handle_check_config(Json(body): Json<Value>) -> Result<Json<Value>, AppError> {
  let request: CheckConfigRequest = serde_json::from_value(body)?;
  Ok(Json(serde_json::to_value(check_config(&request))?))
}

#[utoipa::path(
  post,
  path = "/api/run-config",
  operation_id = "runConfig",
  responses((status = 200, description = "The configuration as a run resolves it, with the outputs the run layer adds and the hash the run records, or the problems found in it"))
)]
async fn handle_run_config(
  State(state): State<Arc<AppState>>,
  Json(body): Json<Value>,
) -> Result<Json<Value>, AppError> {
  let request: RunConfigRequest = serde_json::from_value(body)?;
  let confine = state.confine_hook(request.command)?;
  let response = tokio::task::spawn_blocking(move || run_config(&request, confine)).await?;
  Ok(Json(serde_json::to_value(response)?))
}

#[utoipa::path(
  post,
  path = "/api/check-inputs",
  operation_id = "inputsCheck",
  responses((status = 200, description = "Facts about the input files, read with the readers the commands use"))
)]
async fn handle_check_inputs(
  State(state): State<Arc<AppState>>,
  Json(body): Json<Value>,
) -> Result<Json<Value>, AppError> {
  let mut request: CheckInputsRequest = serde_json::from_value(body)?;
  let policy = state.path_policy()?;
  let confine = |setting: &str, path: &FsPath| {
    policy
      .confine_path(setting, path)
      .map_err(|err| invalid(format!("{err:#}")))
  };
  request.tree = request.tree.as_deref().map(|path| confine("tree", path)).transpose()?;
  request.metadata = request
    .metadata
    .as_deref()
    .map(|path| confine("metadata", path))
    .transpose()?;
  request.alignment = request
    .alignment
    .iter()
    .map(|path| confine("alignment", path))
    .collect::<Result<_, _>>()?;
  let facts = tokio::task::spawn_blocking(move || check_inputs(&request)).await?;
  Ok(Json(serde_json::to_value(facts)?))
}

#[utoipa::path(
  get,
  path = "/api/runs",
  operation_id = "runsList",
  responses((status = 200, description = "Runs, newest first, and the number of runs computing now"))
)]
async fn handle_list_runs(State(state): State<Arc<AppState>>) -> Result<Json<Value>, AppError> {
  Ok(Json(serde_json::to_value(state.runs.list()?)?))
}

#[utoipa::path(
  post,
  path = "/api/runs",
  operation_id = "runsCreate",
  responses((status = 200, description = "The created run; it starts at once unless `defer_start` is set"))
)]
async fn handle_create_run(
  State(state): State<Arc<AppState>>,
  Json(body): Json<Value>,
) -> Result<Json<Value>, AppError> {
  let request: CreateRunRequest = serde_json::from_value(body)?;
  let defer_start = request.defer_start;
  let record = state.runs.create(request)?;
  let record = if defer_start {
    record
  } else {
    start_run(&state, &record.id, None)?
  };
  Ok(Json(serde_json::to_value(record)?))
}

#[utoipa::path(
  get,
  path = "/api/runs/{id}",
  operation_id = "runsGet",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 200, description = "The record of the run"))
)]
async fn handle_get_run(State(state): State<Arc<AppState>>, Path(id): Path<String>) -> Result<Json<Value>, AppError> {
  Ok(Json(serde_json::to_value(state.runs.get(&run_id(&id)?)?)?))
}

#[utoipa::path(
  patch,
  path = "/api/runs/{id}",
  operation_id = "runsUpdate",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 200, description = "The run with its new title or pinned state"))
)]
async fn handle_update_run(
  State(state): State<Arc<AppState>>,
  Path(id): Path<String>,
  Json(body): Json<Value>,
) -> Result<Json<Value>, AppError> {
  let request: UpdateRunRequest = serde_json::from_value(body)?;
  Ok(Json(serde_json::to_value(state.runs.update(&run_id(&id)?, request)?)?))
}

#[utoipa::path(
  delete,
  path = "/api/runs/{id}",
  operation_id = "runsDelete",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 204, description = "The run moved to the trash; restore undoes it"))
)]
async fn handle_delete_run(State(state): State<Arc<AppState>>, Path(id): Path<String>) -> Result<StatusCode, AppError> {
  state.runs.delete(&run_id(&id)?)?;
  Ok(StatusCode::NO_CONTENT)
}

#[utoipa::path(
  post,
  path = "/api/runs/{id}/start",
  operation_id = "runsStart",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 200, description = "The started run"))
)]
async fn handle_start_run(
  State(state): State<Arc<AppState>>,
  Path(id): Path<String>,
  body: Option<Json<Value>>,
) -> Result<Json<Value>, AppError> {
  let request: StartRunRequest = match body {
    Some(Json(body)) => serde_json::from_value(body)?,
    None => StartRunRequest::default(),
  };
  let record = start_run(&state, &run_id(&id)?, request.config)?;
  Ok(Json(serde_json::to_value(record)?))
}

#[utoipa::path(
  post,
  path = "/api/runs/{id}/cancel",
  operation_id = "runsCancel",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 200, description = "Whether cancellation was requested; the run ends with a `cancelled` terminal event", body = CancelRunResponse))
)]
async fn handle_cancel_run(
  State(state): State<Arc<AppState>>,
  Path(id): Path<String>,
) -> Result<Json<CancelRunResponse>, AppError> {
  let cancelled = state.runs.cancel(&run_id(&id)?)?;
  Ok(Json(CancelRunResponse { cancelled }))
}

#[utoipa::path(
  post,
  path = "/api/runs/{id}/restore",
  operation_id = "runsRestore",
  params(("id" = String, Path, description = "Id of a deleted run")),
  responses((status = 200, description = "The run, back from the trash"))
)]
async fn handle_restore_run(
  State(state): State<Arc<AppState>>,
  Path(id): Path<String>,
) -> Result<Json<Value>, AppError> {
  Ok(Json(serde_json::to_value(state.runs.restore(&run_id(&id)?)?)?))
}

#[utoipa::path(
  post,
  path = "/api/runs/{id}/purge",
  operation_id = "runsPurge",
  params(("id" = String, Path, description = "Id of a deleted run")),
  responses((status = 204, description = "The deleted run is removed for good"))
)]
async fn handle_purge_run(State(state): State<Arc<AppState>>, Path(id): Path<String>) -> Result<StatusCode, AppError> {
  state.runs.purge(&run_id(&id)?)?;
  Ok(StatusCode::NO_CONTENT)
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
    (status = 413, description = "The inputs of the run exceed the upload limit of the server", body = ErrorResponse),
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

#[utoipa::path(
  get,
  path = "/api/runs/{id}/files",
  operation_id = "runsFiles",
  params(("id" = String, Path, description = "Id of the run")),
  responses((status = 200, description = "Files in the run's `out/` folder, with their sizes and kinds"))
)]
async fn handle_run_files(State(state): State<Arc<AppState>>, Path(id): Path<String>) -> Result<Json<Value>, AppError> {
  Ok(Json(serde_json::to_value(state.runs.files(&run_id(&id)?)?)?))
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

fn start_run(state: &Arc<AppState>, id: &JobId, config: Option<Value>) -> Result<RunRecord, AppError> {
  let command = state.runs.get(id)?.command;
  let started = state.runs.start(id, config, state.confine_hook(command)?)?;
  let record = started.record().clone();
  let run_id = id.clone();
  drop(tokio::task::spawn_blocking(move || {
    let terminal = started.run();
    match serde_json::to_value(&terminal) {
      Ok(terminal) => info!("Run {} ended: {terminal}", run_id.as_str()),
      Err(err) => error!(
        "Run {} ended; its terminal event cannot be shown: {err}",
        run_id.as_str()
      ),
    }
  }));
  Ok(record)
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

#[derive(OpenApi)]
#[openapi(components(schemas(ErrorResponse, CancelRunResponse)))]
struct SharedSchemas;
