use crate::api::extract::{ApiJson, ApiPath, ApiQuery, OctetStream};
use crate::api::generate::with_project_schemas;
use crate::api::response::{FileContent, TypedSse, ZipAttachment};
use crate::error::AppError;
use crate::events::{app_events_sse, run_events_sse};
use crate::openapi::{add_components, add_discriminators, add_setting_catalog};
use crate::state::{AppState, ServerConfig};
use aide::axum::ApiRouter;
use aide::axum::routing::{get_with, post_with, put_with};
use aide::openapi::{Contact, Info, License, OpenApi};
use app_commands::bridge::service::AppService;
use app_commands::check_config::{CheckConfigRequest, CheckConfigResponse};
use app_commands::check_inputs::{CheckInputsRequest, InputFacts};
use app_commands::datasets::DatasetCatalog;
use app_commands::job::JobId;
use app_commands::results::auspice::AuspiceDocument;
use app_commands::results::clades::{CladeInRuns, CladeRequest};
use app_commands::results::compare::RunComparison;
use app_commands::results::run_results::RunResults;
use app_commands::run_config::{RunConfigRequest, RunConfigResponse};
use app_commands::runs::app_events::AppEvent;
use app_commands::runs::events::RunEvent;
use app_commands::runs::files::RunFile;
use app_commands::runs::manager::UploadedInput;
use app_commands::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use axum::extract::{DefaultBodyLimit, State};
use axum::http::{HeaderMap, HeaderValue};
use axum::response::NoContent;
use axum::{Json, Router};
use eyre::{Report, WrapErr};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::io;
use std::sync::Arc;
use tokio_stream::StreamExt as _;
use tokio_util::io::{StreamReader, SyncIoBridge};
use treetime_schema::{VersionInfo, version_info};

const LAST_EVENT_ID: &str = "last-event-id";

const HEALTH_STATUS: &str = "ok";

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

pub fn api_router(service: Arc<AppService>, config: ServerConfig) -> Result<(Router, OpenApi), Report> {
  let (router, api) = build_api()?;
  let state = Arc::new(AppState {
    config,
    runs: Arc::clone(service.runs()),
    service,
    openapi: serde_json::to_value(&api)?,
  });
  Ok((router.with_state(state), api))
}

pub fn api_doc() -> Result<Value, Report> {
  let (_, api) = build_api()?;
  Ok(serde_json::to_value(api)?)
}

fn build_api() -> Result<(Router<Arc<AppState>>, OpenApi), Report> {
  let mut api = OpenApi {
    info: Info {
      title: "TreeTime API".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
      description: Some(env!("CARGO_PKG_DESCRIPTION").to_owned()),
      contact: Some(Contact {
        name: Some("NeherLab".to_owned()),
        url: Some(env!("CARGO_PKG_HOMEPAGE").to_owned()),
        ..Contact::default()
      }),
      license: Some(License {
        name: env!("CARGO_PKG_LICENSE").to_owned(),
        ..License::default()
      }),
      ..Info::default()
    },
    ..OpenApi::default()
  };
  let router =
    with_project_schemas(|| api_routes().finish_api_with(&mut api, |api| api.default_response::<AppError>()))?;
  add_components(&mut api)?;
  add_discriminators(&mut api)?;
  add_setting_catalog(&mut api)?;
  Ok((router, api))
}

fn api_routes() -> ApiRouter<Arc<AppState>> {
  ApiRouter::new()
    .api_route(
      "/api/health",
      get_with(health, |op| {
        op.id("health")
          .description("Liveness of the server and the version of TreeTime it runs.")
      }),
    )
    .api_route(
      "/api/openapi.json",
      get_with(openapi, |op| op.id("openapi").description("This OpenAPI document.")),
    )
    .api_route(
      "/api/version",
      get_with(version, |op| op.id("version").description("Version of TreeTime.")),
    )
    .api_route(
      "/api/datasets",
      get_with(datasets, |op| {
        op.id("datasets")
          .description("Example datasets and example configurations in the data directory.")
      }),
    )
    .api_route(
      "/api/check-config",
      post_with(config_check, |op| {
        op.id("configCheck")
          .description("The configuration with every default filled in, or the problems found in it.")
      }),
    )
    .api_route(
      "/api/run-config",
      post_with(run_config, |op| {
        op.id("runConfig").description(
          "The configuration as a run resolves it, with the outputs the run layer adds and the hash the run records, \
           or the problems found in it.",
        )
      }),
    )
    .api_route(
      "/api/check-inputs",
      post_with(inputs_check, |op| {
        op.id("inputsCheck")
          .description("Facts about the input files, read with the readers the commands use.")
      }),
    )
    .api_route(
      "/api/runs",
      get_with(runs_list, |op| {
        op.id("runsList")
          .description("Runs, newest first, and the number of runs computing now.")
      })
      .post_with(runs_create, |op| {
        op.id("runsCreate")
          .description("Create a run; it starts at once unless `defer_start` is set.")
      }),
    )
    .api_route(
      "/api/runs/{id}",
      get_with(runs_get, |op| op.id("runsGet").description("The record of a run."))
        .patch_with(runs_update, |op| {
          op.id("runsUpdate")
            .description("Change the title or pinned state of a run.")
        })
        .delete_with(runs_delete, |op| {
          op.id("runsDelete")
            .description("Move a run to the trash; restore undoes it.")
        }),
    )
    .api_route(
      "/api/runs/{id}/start",
      post_with(runs_start, |op| op.id("runsStart").description("Start a created run.")),
    )
    .api_route(
      "/api/runs/{id}/cancel",
      post_with(runs_cancel, |op| {
        op.id("runsCancel")
          .description("Request cancellation of a run; the run ends with a `cancelled` terminal event.")
      }),
    )
    .api_route(
      "/api/runs/{id}/restore",
      post_with(runs_restore, |op| {
        op.id("runsRestore").description("Bring a run back from the trash.")
      }),
    )
    .api_route(
      "/api/runs/{id}/purge",
      post_with(runs_purge, |op| {
        op.id("runsPurge").description("Remove a deleted run for good.")
      }),
    )
    .api_route(
      "/api/runs/{id}/files",
      get_with(runs_files, |op| {
        op.id("runsFiles")
          .description("Files in the run's `out/` folder, with their sizes and kinds.")
      }),
    )
    .api_route(
      "/api/runs/{id}/results",
      get_with(runs_results, |op| {
        op.id("runsResults")
          .description("Results of a finished run, read from its output files.")
      }),
    )
    .api_route(
      "/api/runs/{id}/auspice",
      get_with(runs_auspice, |op| {
        op.id("runsAuspice")
          .description("Auspice JSON of a finished run, with the color scales the app displays.")
      }),
    )
    .api_route(
      "/api/runs/{id}/compare/{other}",
      get_with(runs_compare, |op| {
        op.id("runsCompare")
          .description("Differences of the second run's results from the first run's.")
      }),
    )
    .api_route(
      "/api/clade-in-runs",
      post_with(clade_in_runs, |op| {
        op.id("cladeInRuns")
          .description("Nodes of the other finished time-tree runs with the same set of samples below them.")
      }),
    )
    .api_route(
      "/api/runs/{id}/events",
      get_with(runs_events, |op| {
        op.id("runsEvents").description(
          "Stream of the run's events from `from`: `started`, `progress`, `log` and `iteration`, then one \
           `terminal`. The `Last-Event-ID` header of a reconnect resumes after that event.",
        )
      }),
    )
    .api_route(
      "/api/events",
      get_with(events, |op| {
        op.id("events").description(
          "Stream of changes to the runs: `run-created`, `run-updated`, `run-deleted`, `run-restored` and \
           `run-purged`, each with the REST paths it made stale, either the path alone (`exact`) or the path and \
           every path below it (`subtree`). Without `from` the stream sends the changes from now \
           on. `from`, or the `Last-Event-ID` header of a reconnect, resumes after an earlier event; when the server \
           no longer keeps that event or the event is from a previous server, the stream starts with a `resync` \
           event instead.",
        )
      }),
    )
    .api_route(
      "/api/runs/{id}/inputs/{name}",
      put_with(runs_upload_input, |op| {
        op.id("runsUploadInput")
          .description(
            "Store a file in the `inputs/` folder of a run that has not started; the answer names the path to use \
             for it in the run's configuration.",
          )
          .response_with::<413, AppError, _>(|response| {
            response.description("The inputs of the run exceed the upload limit of the server")
          })
      })
      .layer(DefaultBodyLimit::disable()),
    )
    .api_route(
      "/api/runs/{id}/file",
      get_with(runs_file, |op| {
        op.id("runsFile")
          .description("Contents of a file in the run's `out/` folder.")
      }),
    )
    .api_route(
      "/api/runs/{id}/archive",
      get_with(runs_archive, |op| {
        op.id("runsArchive")
          .description("Zip archive of the run's `out/` folder.")
      }),
    )
}

/// Liveness of the server.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
struct HealthStatus {
  /// Always `ok`.
  status: String,
  /// Version of TreeTime.
  version: String,
}

/// Path of a run.
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[schemars(inline)]
struct RunPath {
  /// Id of the run.
  id: JobId,
}

/// Path of two runs.
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[schemars(inline)]
struct RunPairPath {
  /// Id of the first run.
  id: JobId,
  /// Id of the run compared with the first.
  other: JobId,
}

/// Path of an input file of a run.
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[schemars(inline)]
struct InputPath {
  /// Id of a run that has not started.
  id: JobId,
  /// File name inside the run's `inputs/` folder.
  name: String,
}

/// Query of the run event stream.
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[schemars(inline)]
struct EventsQuery {
  /// Sequence number of the first event to send.
  from: Option<usize>,
}

/// Query of a run file.
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[schemars(inline)]
struct FileQuery {
  /// Path of the file relative to the run's `out/` folder.
  path: String,
}

async fn health() -> Json<HealthStatus> {
  Json(HealthStatus {
    status: HEALTH_STATUS.to_owned(),
    version: version_info().version.to_owned(),
  })
}

async fn openapi(State(state): State<Arc<AppState>>) -> Json<Value> {
  Json(state.openapi.clone())
}

async fn version(State(state): State<Arc<AppState>>) -> Result<Json<VersionInfo>, AppError> {
  call(&state, AppService::version).await.map(Json)
}

async fn datasets(State(state): State<Arc<AppState>>) -> Result<Json<DatasetCatalog>, AppError> {
  call(&state, AppService::datasets).await.map(Json)
}

async fn config_check(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CheckConfigRequest>,
) -> Result<Json<CheckConfigResponse>, AppError> {
  call(&state, move |service| service.check_config(&request))
    .await
    .map(Json)
}

async fn run_config(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<RunConfigRequest>,
) -> Result<Json<RunConfigResponse>, AppError> {
  call(&state, move |service| service.run_config(&request))
    .await
    .map(Json)
}

async fn inputs_check(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CheckInputsRequest>,
) -> Result<Json<InputFacts>, AppError> {
  call(&state, move |service| service.check_inputs(request))
    .await
    .map(Json)
}

async fn runs_list(State(state): State<Arc<AppState>>) -> Result<Json<RunList>, AppError> {
  call(&state, AppService::list_runs).await.map(Json)
}

async fn runs_create(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CreateRunRequest>,
) -> Result<Json<RunRecord>, AppError> {
  call(&state, move |service| service.create_run(request)).await.map(Json)
}

async fn runs_get(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<Json<RunRecord>, AppError> {
  call(&state, move |service| service.get_run(&id)).await.map(Json)
}

async fn runs_start(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiJson(request): ApiJson<StartRunRequest>,
) -> Result<Json<RunRecord>, AppError> {
  call(&state, move |service| service.start_run(&id, request))
    .await
    .map(Json)
}

async fn runs_update(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiJson(request): ApiJson<UpdateRunRequest>,
) -> Result<Json<RunSummary>, AppError> {
  call(&state, move |service| service.update_run(&id, request))
    .await
    .map(Json)
}

async fn runs_cancel(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<Json<CancelRunResponse>, AppError> {
  call(&state, move |service| service.cancel_run(&id)).await.map(Json)
}

async fn runs_delete(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<NoContent, AppError> {
  call(&state, move |service| service.delete_run(&id)).await?;
  Ok(NoContent)
}

async fn runs_restore(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<Json<RunSummary>, AppError> {
  call(&state, move |service| service.restore_run(&id)).await.map(Json)
}

async fn runs_purge(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<NoContent, AppError> {
  call(&state, move |service| service.purge_run(&id)).await?;
  Ok(NoContent)
}

async fn runs_files(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<Json<Vec<RunFile>>, AppError> {
  call(&state, move |service| service.run_files(&id)).await.map(Json)
}

async fn runs_results(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<Json<RunResults>, AppError> {
  call(&state, move |service| service.run_results(&id)).await.map(Json)
}

async fn runs_auspice(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<Json<AuspiceDocument>, AppError> {
  call(&state, move |service| service.run_auspice(&id)).await.map(Json)
}

async fn runs_compare(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPairPath { id, other }): ApiPath<RunPairPath>,
) -> Result<Json<RunComparison>, AppError> {
  call(&state, move |service| service.compare_runs(&id, &other))
    .await
    .map(Json)
}

async fn clade_in_runs(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CladeRequest>,
) -> Result<Json<CladeInRuns>, AppError> {
  call(&state, move |service| service.clade_in_runs(&request))
    .await
    .map(Json)
}

async fn runs_events(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiQuery(query): ApiQuery<EventsQuery>,
  headers: HeaderMap,
) -> Result<TypedSse<RunEvent>, AppError> {
  run_events_sse(&state, &id, resume_from(&headers, &query).unwrap_or(0))
}

async fn events(
  State(state): State<Arc<AppState>>,
  ApiQuery(query): ApiQuery<EventsQuery>,
  headers: HeaderMap,
) -> TypedSse<AppEvent> {
  app_events_sse(&state, resume_from(&headers, &query))
}

fn resume_from(headers: &HeaderMap, query: &EventsQuery) -> Option<usize> {
  headers
    .get(LAST_EVENT_ID)
    .and_then(|value| value.to_str().ok())
    .and_then(|value| value.parse::<usize>().ok())
    .map(|seq| seq + 1)
    .or(query.from)
}

async fn runs_upload_input(
  State(state): State<Arc<AppState>>,
  ApiPath(InputPath { id, name }): ApiPath<InputPath>,
  OctetStream(body): OctetStream,
) -> Result<Json<UploadedInput>, AppError> {
  let stream = body.into_data_stream().map(|chunk| chunk.map_err(io::Error::other));
  let mut reader = SyncIoBridge::new(StreamReader::new(stream));
  let runs = Arc::clone(&state.runs);
  let limit = state.config.max_upload_size;
  let uploaded = tokio::task::spawn_blocking(move || runs.upload_input(&id, &name, &mut reader, limit)).await??;
  Ok(Json(uploaded))
}

async fn runs_file(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiQuery(query): ApiQuery<FileQuery>,
) -> Result<FileContent, AppError> {
  let path = state.runs.file_path(&id, &query.path)?;
  let bytes = tokio::fs::read(&path)
    .await
    .wrap_err_with(|| format!("When reading '{}'", path.display()))?;
  let extension = path.extension().and_then(|extension| extension.to_str()).unwrap_or("");
  Ok(FileContent {
    content_type: content_type(extension),
    bytes,
  })
}

async fn runs_archive(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ZipAttachment, AppError> {
  let runs = Arc::clone(&state.runs);
  let archive_id = id.clone();
  let bytes = tokio::task::spawn_blocking(move || runs.zip(&archive_id)).await??;
  let disposition = HeaderValue::from_str(&format!("attachment; filename=\"treetime-{}.zip\"", id.as_str()))?;
  Ok(ZipAttachment { disposition, bytes })
}

async fn call<T, F>(state: &AppState, operation: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce(&AppService) -> Result<T, Report> + Send + 'static,
{
  let service = Arc::clone(&state.service);
  Ok(tokio::task::spawn_blocking(move || operation(&service)).await??)
}

fn content_type(extension: &str) -> HeaderValue {
  let content_type = CONTENT_TYPES
    .iter()
    .find(|(known, _)| *known == extension)
    .map_or("application/octet-stream", |(_, content_type)| *content_type);
  HeaderValue::from_static(content_type)
}
