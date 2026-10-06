use crate::api::cache::{Revalidated, finished_runs_tag, is_unchanged};
use crate::api::download::{serve_run_file, stream_run_archive};
use crate::api::extract::{ApiJson, ApiPath, ApiQuery, OctetStream};
use crate::api::generate::with_project_schemas;
use crate::api::response::{FileContent, TypedSse, ZipAttachment};
use crate::app_settings_routes::{app_settings_routes, host_routes};
use crate::error::{AppError, panic_response, plain_error};
use crate::events::{app_events_sse, run_events_sse};
use crate::openapi::{add_components, add_discriminators};
use crate::state::{AppState, ServerConfig};
use aide::axum::ApiRouter;
use aide::axum::routing::{get_with, post_with, put_with};
use aide::openapi::{Contact, Info, License, OpenApi};
use app_commands::bridge::error::ErrorCode;
use app_commands::bridge::service::AppService;
use app_commands::check_config::{CheckConfigRequest, CheckConfigResponse};
use app_commands::check_inputs::{CheckInputsRequest, InputFacts};
use app_commands::config::schema::draft2020_settings;
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
use axum::error_handling::HandleErrorLayer;
use axum::extract::{DefaultBodyLimit, Request, State};
use axum::http::{HeaderMap, HeaderValue, Method, Uri, header};
use axum::response::{IntoResponse, Response};
use axum::{BoxError, Router};
use eyre::{Report, eyre};
use itertools::Itertools;
use schemars::JsonSchema;
use schemars::generate::SchemaSettings;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::any::Any;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io;
use std::process;
use std::sync::Arc;
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use tokio_stream::StreamExt as _;
use tokio_util::io::{StreamReader, SyncIoBridge};
use tower::ServiceBuilder;
use tower::timeout::error::Elapsed;
use tower_http::catch_panic::CatchPanicLayer;
use tower_http::set_header::SetResponseHeaderLayer;
use treetime_schema::{VersionInfo, version_info};

const LAST_EVENT_ID: &str = "last-event-id";

const HEALTH_STATUS: &str = "ok";

pub(crate) const NO_STORE: &str = "no-store";

const REQUEST_TIMEOUT: Duration = Duration::from_secs(120);

pub fn api_router(service: Arc<AppService>, config: ServerConfig) -> Result<(Router, OpenApi), Report> {
  let scope = if config.settings.is_some() {
    RouteScope::Renderer
  } else {
    RouteScope::Public
  };
  let (routes, api) = build_api(scope, &draft2020_settings())?;
  let state = app_state(service, config, &api)?;
  Ok((finish_router(routes, state), api))
}

pub fn local_api_routers(service: Arc<AppService>, config: ServerConfig) -> Result<LocalRouters, Report> {
  let (host, api) = build_api(RouteScope::Host, &draft2020_settings())?;
  let (renderer, _) = build_api(RouteScope::Renderer, &draft2020_settings())?;
  let state = app_state(service, config, &api)?;
  Ok(LocalRouters {
    host: finish_router(host, Arc::clone(&state)),
    renderer: finish_router(renderer, state),
  })
}

pub fn api_doc() -> Result<Value, Report> {
  api_doc_with(&draft2020_settings())
}

pub(crate) fn api_doc_with(settings: &SchemaSettings) -> Result<Value, Report> {
  let (_, api) = build_api(RouteScope::Host, settings)?;
  Ok(serde_json::to_value(api)?)
}

/// Routers of a local app: the host router also saves run files to a path that the caller names, so only the process
/// that hosts the app may reach it.
pub struct LocalRouters {
  pub host: Router,
  pub renderer: Router,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RouteScope {
  Public,
  Renderer,
  Host,
}

fn app_state(service: Arc<AppService>, config: ServerConfig, api: &OpenApi) -> Result<Arc<AppState>, Report> {
  Ok(Arc::new(AppState {
    config,
    runs: Arc::clone(service.runs()),
    service,
    openapi: serde_json::to_value(api)?,
    instance: instance_id()?,
  }))
}

fn finish_router(routes: Router<Arc<AppState>>, state: Arc<AppState>) -> Router {
  routes
    .method_not_allowed_fallback(method_not_allowed)
    .fallback(not_found)
    .layer(SetResponseHeaderLayer::if_not_present(
      header::CACHE_CONTROL,
      HeaderValue::from_static(NO_STORE),
    ))
    .layer(CatchPanicLayer::custom(|payload: Box<dyn Any + Send>| {
      panic_response(&*payload)
    }))
    .with_state(state)
}

fn build_api(scope: RouteScope, settings: &SchemaSettings) -> Result<(Router<Arc<AppState>>, OpenApi), Report> {
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
  let router = with_project_schemas(settings, || {
    let routes = match scope {
      RouteScope::Public => api_routes(),
      RouteScope::Renderer => api_routes().merge(app_settings_routes()),
      RouteScope::Host => api_routes().merge(app_settings_routes()).merge(host_routes()),
    };
    routes.finish_api_with(&mut api, |api| api.default_response::<AppError>())
  })?;
  add_components(&mut api, settings)?;
  add_discriminators(&mut api)?;
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
          .description("Example datasets and example configurations in the examples folder.")
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
          "Stream of changes to the runs: `run-created` and `run-updated`, each with the REST paths it made stale, \
           either the path alone (`exact`) or the path and every path below it (`subtree`). Without `from` the stream sends the changes from now \
           on. `from`, or the `Last-Event-ID` header of a reconnect, resumes after an earlier event; when the server \
           no longer keeps that event or the event is from a previous server, the stream starts with a `resync` \
           event instead.",
        )
      }),
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
    .layer(
      ServiceBuilder::new()
        .layer(HandleErrorLayer::new(timed_out))
        .timeout(REQUEST_TIMEOUT),
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
}

/// Liveness of the server.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
struct HealthStatus {
  /// Always `ok`.
  status: String,
  /// Version of TreeTime.
  version: String,
}

/// Path of a run.
#[derive(Clone, Debug, Deserialize, JsonSchema, deser::Deserialize)]
#[schemars(inline)]
pub(crate) struct RunPath {
  /// Id of the run.
  pub(crate) id: JobId,
}

/// Path of two runs.
#[derive(Clone, Debug, Deserialize, JsonSchema, deser::Deserialize)]
#[schemars(inline)]
struct RunPairPath {
  /// Id of the first run.
  id: JobId,
  /// Id of the run compared with the first.
  other: JobId,
}

/// Path of an input file of a run.
#[derive(Clone, Debug, Deserialize, JsonSchema, deser::Deserialize)]
#[schemars(inline)]
struct InputPath {
  /// Id of a run that has not started.
  id: JobId,
  /// File name inside the run's `inputs/` folder.
  name: String,
}

/// Query of the run event stream.
#[derive(Clone, Debug, Deserialize, JsonSchema, deser::Deserialize)]
#[schemars(inline)]
struct EventsQuery {
  /// Sequence number of the first event to send.
  from: Option<usize>,
}

/// Query of a run file.
#[derive(Clone, Debug, Deserialize, JsonSchema, deser::Deserialize)]
#[schemars(inline)]
struct FileQuery {
  /// Path of the file relative to the run's `out/` folder.
  path: String,
}

async fn health() -> ApiJson<HealthStatus> {
  ApiJson(HealthStatus {
    status: HEALTH_STATUS.to_owned(),
    version: version_info().version.to_owned(),
  })
}

async fn openapi(State(state): State<Arc<AppState>>) -> axum::Json<Value> {
  axum::Json(state.openapi.clone())
}

async fn version(State(state): State<Arc<AppState>>) -> Result<ApiJson<VersionInfo>, AppError> {
  call(&state, AppService::version).await.map(ApiJson)
}

async fn datasets(State(state): State<Arc<AppState>>) -> Result<ApiJson<DatasetCatalog>, AppError> {
  call(&state, AppService::datasets).await.map(ApiJson)
}

async fn config_check(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CheckConfigRequest>,
) -> Result<ApiJson<CheckConfigResponse>, AppError> {
  call(&state, move |service| service.check_config(&request))
    .await
    .map(ApiJson)
}

async fn run_config(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<RunConfigRequest>,
) -> Result<ApiJson<RunConfigResponse>, AppError> {
  call(&state, move |service| service.run_config(&request))
    .await
    .map(ApiJson)
}

async fn inputs_check(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CheckInputsRequest>,
) -> Result<ApiJson<InputFacts>, AppError> {
  call(&state, move |service| service.check_inputs(request))
    .await
    .map(ApiJson)
}

async fn runs_list(State(state): State<Arc<AppState>>) -> Result<ApiJson<RunList>, AppError> {
  call(&state, AppService::list_runs).await.map(ApiJson)
}

async fn runs_create(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CreateRunRequest>,
) -> Result<ApiJson<RunRecord>, AppError> {
  call(&state, move |service| service.create_run(&request))
    .await
    .map(ApiJson)
}

async fn runs_get(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ApiJson<RunRecord>, AppError> {
  call(&state, move |service| service.get_run(&id)).await.map(ApiJson)
}

async fn runs_start(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiJson(request): ApiJson<StartRunRequest>,
) -> Result<ApiJson<RunRecord>, AppError> {
  call(&state, move |service| service.start_run(&id, &request))
    .await
    .map(ApiJson)
}

async fn runs_update(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiJson(request): ApiJson<UpdateRunRequest>,
) -> Result<ApiJson<RunSummary>, AppError> {
  call(&state, move |service| service.update_run(&id, request))
    .await
    .map(ApiJson)
}

async fn runs_cancel(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ApiJson<CancelRunResponse>, AppError> {
  call(&state, move |service| service.cancel_run(&id)).await.map(ApiJson)
}

async fn runs_files(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ApiJson<Vec<RunFile>>, AppError> {
  call(&state, move |service| service.run_files(&id)).await.map(ApiJson)
}

async fn runs_results(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  headers: HeaderMap,
) -> Result<Revalidated<ApiJson<RunResults>>, AppError> {
  revalidated(&state, &headers, vec![id.clone()], move |service| {
    service.run_results(&id)
  })
  .await
}

async fn runs_auspice(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  headers: HeaderMap,
) -> Result<Revalidated<ApiJson<AuspiceDocument>>, AppError> {
  revalidated(&state, &headers, vec![id.clone()], move |service| {
    service.run_auspice(&id)
  })
  .await
}

async fn runs_compare(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPairPath { id, other }): ApiPath<RunPairPath>,
  headers: HeaderMap,
) -> Result<Revalidated<ApiJson<RunComparison>>, AppError> {
  revalidated(&state, &headers, vec![id.clone(), other.clone()], move |service| {
    service.compare_runs(&id, &other)
  })
  .await
}

async fn clade_in_runs(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<CladeRequest>,
) -> Result<ApiJson<CladeInRuns>, AppError> {
  call(&state, move |service| service.clade_in_runs(&request))
    .await
    .map(ApiJson)
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
) -> Result<ApiJson<UploadedInput>, AppError> {
  let stream = body.into_data_stream().map(|chunk| chunk.map_err(io::Error::other));
  let mut reader = SyncIoBridge::new(StreamReader::new(stream));
  let runs = Arc::clone(&state.runs);
  let limit = state.config.max_upload_size;
  let uploaded = tokio::task::spawn_blocking(move || runs.upload_input(&id, &name, &mut reader, limit)).await??;
  Ok(ApiJson(uploaded))
}

async fn runs_file(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiQuery(query): ApiQuery<FileQuery>,
  request: Request,
) -> Result<FileContent, AppError> {
  let path = call(&state, move |service| service.runs().file_path(&id, &query.path)).await?;
  serve_run_file(path, request).await
}

async fn runs_archive(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ZipAttachment, AppError> {
  let name = id.as_str().to_owned();
  let out_dir = call(&state, move |service| service.runs().out_dir(&id)).await?;
  stream_run_archive(out_dir, &name)
}

async fn revalidated<T, F>(
  state: &AppState,
  headers: &HeaderMap,
  ids: Vec<JobId>,
  operation: F,
) -> Result<Revalidated<ApiJson<T>>, AppError>
where
  T: Send + 'static,
  F: FnOnce(&AppService) -> Result<T, Report> + Send + 'static,
{
  let instance = state.instance;
  let tag = call(state, move |service| {
    let records = ids
      .iter()
      .map(|id| service.runs().get(id))
      .try_collect::<_, Vec<_>, _>()?;
    finished_runs_tag(instance, &records)
  })
  .await?;
  if let Some(tag) = tag.clone().filter(|tag| is_unchanged(headers, Some(tag))) {
    return Ok(Revalidated::Unchanged(tag));
  }
  let value = call(state, operation).await?;
  Ok(Revalidated::Fresh(tag, ApiJson(value)))
}

async fn call<T, F>(state: &AppState, operation: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce(&AppService) -> Result<T, Report> + Send + 'static,
{
  let service = Arc::clone(&state.service);
  Ok(tokio::task::spawn_blocking(move || operation(&service)).await??)
}

pub(crate) async fn not_found(uri: Uri) -> Response {
  plain_error(ErrorCode::NotFound, format!("the API has no path `{}`", uri.path()))
}

async fn method_not_allowed(method: Method, uri: Uri) -> Response {
  plain_error(
    ErrorCode::MethodNotAllowed,
    format!("the API path `{}` does not accept `{method}` requests", uri.path()),
  )
}

async fn timed_out(error: BoxError) -> Response {
  if error.is::<Elapsed>() {
    plain_error(
      ErrorCode::Timeout,
      format!("the request took longer than {} seconds", REQUEST_TIMEOUT.as_secs()),
    )
  } else {
    AppError::from(eyre!(error)).into_response()
  }
}

fn instance_id() -> Result<u64, Report> {
  let mut hasher = DefaultHasher::new();
  SystemTime::now()
    .duration_since(UNIX_EPOCH)?
    .as_nanos()
    .hash(&mut hasher);
  process::id().hash(&mut hasher);
  Ok(hasher.finish())
}
