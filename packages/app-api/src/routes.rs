use crate::api::blocking::{receive_body, run_blocking};
use crate::api::download::{serve_run_file, stream_run_package};
use crate::api::extract::{ApiJson, ApiPath, ApiQuery, OctetStream};
use crate::api::generate::with_project_schemas;
use crate::api::response::{FileContent, TypedSse, ZipAttachment};
use crate::app_settings_routes::{app_settings_routes, host_routes};
use crate::error::{AppError, panic_response, plain_error};
use crate::events::{app_events_sse, run_events_sse};
use crate::openapi::{add_components, add_discriminators, remove_cli_only_commands};
use crate::state::{AppState, LocalSettings, RouterConfig};
use aide::axum::ApiRouter;
use aide::axum::routing::{get_with, post_with, put_with};
use aide::openapi::{Contact, Info, License, OpenApi};
use app_backend::bridge::error::ErrorCode;
use app_backend::bridge::service::AppService;
use app_backend::check_config::{CheckConfigRequest, CheckConfigResponse};
use app_backend::check_inputs::{CheckInputsRequest, InputFacts};
use app_backend::clock_preview::{ClockPreview, ClockPreviewRequest};
use app_backend::file_kinds::{FileKinds, FileKindsRequest};
use app_backend::job::JobId;
use app_backend::results::auspice::AuspiceDocument;
use app_backend::results::clades::{CladeInRuns, CladeRequest};
use app_backend::results::compare::RunComparison;
use app_backend::results::run_results::RunResults;
use app_backend::run_config::{RunConfigRequest, RunConfigResponse};
use app_backend::runs::app_events::AppEvent;
use app_backend::runs::events::RunEvent;
use app_backend::runs::files::RunFile;
use app_backend::runs::manager::{UploadAllowance, UploadedInput};
use app_backend::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use app_commands::config::schema::draft2020_settings;
use app_commands::json_value::SparseConfig;
use axum::Router;
use axum::extract::{DefaultBodyLimit, Request, State};
use axum::http::{HeaderMap, HeaderValue, Method, Uri, header};
use axum::response::Response;
use deser::Deserialize;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use schemars::generate::SchemaSettings;
use serde_json::Value;
use std::any::Any;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tower_http::catch_panic::CatchPanicLayer;
use tower_http::set_header::SetResponseHeaderLayer;
use treetime_schema::VersionInfo;

const LAST_EVENT_ID: &str = "last-event-id";

pub(crate) const NO_STORE: &str = "no-store";

pub fn api_router(service: Arc<AppService>, config: RouterConfig) -> Result<(Router, OpenApi), Report> {
  let (routes, api) = build_api(RouteScope::Public, &draft2020_settings())?;
  let state = app_state(service, config, None);
  Ok((finish_router(routes, state), api))
}

pub fn local_api_routers(
  service: Arc<AppService>,
  config: RouterConfig,
  settings: LocalSettings,
) -> Result<LocalRouters, Report> {
  let (host, _) = build_api(RouteScope::Host, &draft2020_settings())?;
  let (renderer, _) = build_api(RouteScope::Renderer, &draft2020_settings())?;
  let state = app_state(service, config, Some(settings));
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

pub struct LocalRouters {
  pub host: Router,
  pub renderer: Router,
}

fn app_state(service: Arc<AppService>, config: RouterConfig, settings: Option<LocalSettings>) -> Arc<AppState> {
  Arc::new(AppState {
    config,
    settings,
    runs: Arc::clone(service.runs()),
    service,
  })
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
  remove_cli_only_commands(&mut api)?;
  add_discriminators(&mut api)?;
  Ok((router, api))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RouteScope {
  Public,
  Renderer,
  Host,
}

fn api_routes() -> ApiRouter<Arc<AppState>> {
  ApiRouter::new()
    .api_route(
      "/api/version",
      get_with(version, |op| op.id("version").description("Version of TreeTime.")),
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
      "/api/file-kinds",
      post_with(file_kinds, |op| {
        op.id("fileKinds").description(
          "The kind of each file (tree, alignment, metadata, configuration, run package), read from its content with \
           the readers the commands use, or `unknown` with the reason.",
        )
      }),
    )
    .api_route(
      "/api/clock-preview",
      post_with(clock_preview, |op| {
        op.id("clockPreview").description(
          "The root-to-tip regression of a draft's inputs, rooting and clock filter, computed as a `clock` run \
           computes it, without a run, and the tree in divergence units; trees above the tip limit only on request.",
        )
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
          op.id("runsDelete").description(
            "Delete a run that is not computing, with its folder and the uploaded inputs that no other run \
             references. Cancel a computing run first. A run that is already deleted answers not found.",
          )
        }),
    )
    .api_route(
      "/api/runs/{id}/settings",
      get_with(runs_settings, |op| {
        op.id("runsSettings").description(
          "Settings of a run, with each input path as requests take it and without the output paths, to start a \
           new run from the settings of this one.",
        )
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
      "/api/runs/{id}/inputs",
      get_with(runs_upload_allowance, |op| {
        op.id("runsUploadAllowance").description(
          "Report how many bytes the next uploads of a run that has not started may add before the inputs of the run \
           reach the upload limit.",
        )
      })
      .delete_with(runs_release_input, |op| {
        op.id("runsReleaseInput").description(
          "Remove an uploaded file from a run that has not started; the file stays while another run uses it. \
           Removing an upload that is already gone succeeds.",
        )
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
      "/api/compare/{id}/{other}",
      get_with(runs_compare, |op| {
        op.id("runsCompare").description(
          "Differences of the second run's settings and results from the first run's. Every change of a run makes \
           the comparisons under `/api/compare` stale.",
        )
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
          "Stream of changes to the runs: `run-created`, `run-updated`, and `run-removed`, each with the REST paths \
           it made stale, either the path alone (`exact`) or the path and every path below it (`subtree`). Without \
           `from` the stream sends the changes from now on. `from`, or the `Last-Event-ID` header of a reconnect, resumes after an earlier event; when the back end \
           no longer keeps that event or the event is from an earlier back end, the stream starts with a `resync` \
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
      "/api/runs/{id}/package",
      get_with(runs_package, |op| {
        op.id("runsPackage").description(
          "Run package of a run that has ended: a zip archive with one folder `treetime-<id>/` that holds `run.json`, \
           `events.jsonl`, the input files in `inputs/`, the output files in `out/`, the YAML config `run.yaml`, and \
           `command.sh`, which runs `run.yaml` with the command-line tool. Every path in `run.json` and `run.yaml` is \
           relative to the folder.",
        )
      }),
    )
    .api_route(
      "/api/runs/import",
      post_with(runs_import, |op| {
        op.id("runsImport").description(
          "Create a run from a run package. The package must hold a run that has ended, every path in it must be \
           relative and inside it, every input must match the size and SHA-256 that `run.json` lists, and every \
           output file that `run.json` lists must be present. The run keeps the id of the package unless a run with \
           that id exists. The package, as uploaded and unpacked, must fit in the upload limit.",
        )
        .response_with::<413, AppError, _>(|response| {
          response.description("The run package exceeds the upload limit")
        })
      })
      .layer(DefaultBodyLimit::disable()),
    )
    .api_route(
      "/api/runs/{id}/inputs/{name}",
      put_with(runs_upload_input, |op| {
        op.id("runsUploadInput")
          .description(
            "Store a file for a run that has not started; the answer names the path to use for it in the run's \
             configuration.",
          )
          .response_with::<413, AppError, _>(|response| {
            response.description("The inputs of the run exceed the upload limit")
          })
      })
      .layer(DefaultBodyLimit::disable()),
    )
}

async fn version(State(state): State<Arc<AppState>>) -> Result<ApiJson<VersionInfo>, AppError> {
  call(&state, AppService::version).await.map(ApiJson)
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

async fn file_kinds(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<FileKindsRequest>,
) -> Result<ApiJson<FileKinds>, AppError> {
  call(&state, move |service| service.file_kinds(&request))
    .await
    .map(ApiJson)
}

async fn clock_preview(
  State(state): State<Arc<AppState>>,
  ApiJson(request): ApiJson<ClockPreviewRequest>,
) -> Result<ApiJson<ClockPreview>, AppError> {
  call(&state, move |service| service.clock_preview(request))
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

async fn runs_delete(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<(), AppError> {
  call(&state, move |service| service.delete_run(&id)).await
}

async fn runs_settings(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ApiJson<SparseConfig>, AppError> {
  call(&state, move |service| service.run_settings(&id))
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
) -> Result<ApiJson<RunResults>, AppError> {
  call(&state, move |service| service.run_results(&id)).await.map(ApiJson)
}

async fn runs_auspice(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ApiJson<AuspiceDocument>, AppError> {
  call(&state, move |service| service.run_auspice(&id)).await.map(ApiJson)
}

async fn runs_compare(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPairPath { id, other }): ApiPath<RunPairPath>,
) -> Result<ApiJson<RunComparison>, AppError> {
  call(&state, move |service| service.compare_runs(&id, &other))
    .await
    .map(ApiJson)
}

/// Path of two runs.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
struct RunPairPath {
  /// Id of the first run.
  id: JobId,
  /// Id of the run compared with the first.
  other: JobId,
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

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    result_defaulted,
    reason = "a `Last-Event-ID` header that does not parse as an event number gives no resume point, as an absent one"
  )
)]
fn resume_from(headers: &HeaderMap, query: &EventsQuery) -> Option<u64> {
  headers
    .get(LAST_EVENT_ID)
    .and_then(|value| value.to_str().ok())
    .and_then(|value| value.parse::<u64>().ok())
    .map(|seq| seq + 1)
    .or(query.from)
}

/// Query of the run event stream.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
struct EventsQuery {
  /// Sequence number of the first event to send.
  from: Option<u64>,
}

async fn runs_upload_input(
  State(state): State<Arc<AppState>>,
  ApiPath(InputPath { id, name }): ApiPath<InputPath>,
  ApiQuery(UploadInputQuery { replaces }): ApiQuery<UploadInputQuery>,
  OctetStream(body): OctetStream,
) -> Result<ApiJson<UploadedInput>, AppError> {
  let runs = Arc::clone(&state.runs);
  let limit = state.config.max_upload_size;
  let replaces = replaces.into_iter().map(PathBuf::from).collect_vec();
  let upload = run_blocking(move || runs.begin_upload(&id, &name, &replaces, limit)).await?;
  Ok(ApiJson(receive_body(body, upload).await?))
}

async fn runs_upload_allowance(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiQuery(UploadInputQuery { replaces }): ApiQuery<UploadInputQuery>,
) -> Result<ApiJson<UploadAllowance>, AppError> {
  let runs = Arc::clone(&state.runs);
  let limit = state.config.max_upload_size;
  let replaces = replaces.into_iter().map(PathBuf::from).collect_vec();
  run_blocking(move || runs.upload_allowance(&id, &replaces, limit))
    .await
    .map(ApiJson)
}

/// Path of an input file of a run.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
struct InputPath {
  /// Id of a run that has not started.
  id: JobId,
  /// File name of the upload.
  name: String,
}

/// Query of an upload or of the upload allowance of a run.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
struct UploadInputQuery {
  /// Paths of earlier uploads of the run that the upload replaces, as their upload answers gave them; their size
  /// does not count toward the upload limit.
  #[deser(default)]
  #[schemars(default)]
  replaces: Vec<String>,
}

async fn runs_release_input(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiQuery(query): ApiQuery<ReleaseInputQuery>,
) -> Result<(), AppError> {
  call(&state, move |service| {
    service.runs().release_input(&id, Path::new(&query.path))
  })
  .await
}

/// Query of an upload to remove.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
struct ReleaseInputQuery {
  /// Path of the upload, as the upload answer gave it.
  path: String,
}

async fn runs_file(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiQuery(query): ApiQuery<FileQuery>,
  request: Request,
) -> Result<FileContent, AppError> {
  let download = call(&state, move |service| service.runs().file_download(&id, &query.path)).await?;
  serve_run_file(download, request).await
}

/// Query of a run file.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
struct FileQuery {
  /// Path of the file relative to the run's `out/` folder.
  path: String,
}

async fn runs_package(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
) -> Result<ZipAttachment, AppError> {
  let package = call(&state, move |service| service.runs().package(&id)).await?;
  stream_run_package(package)
}

async fn runs_import(
  State(state): State<Arc<AppState>>,
  OctetStream(body): OctetStream,
) -> Result<ApiJson<RunRecord>, AppError> {
  let runs = Arc::clone(&state.runs);
  let limit = state.config.max_upload_size;
  let import = run_blocking(move || runs.begin_import(limit)).await?;
  Ok(ApiJson(receive_body(body, import).await?))
}

/// Path of a run.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(inline)]
pub(crate) struct RunPath {
  /// Id of the run.
  pub(crate) id: JobId,
}

async fn call<T, F>(state: &AppState, operation: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce(&AppService) -> Result<T, Report> + Send + 'static,
{
  let service = Arc::clone(&state.service);
  run_blocking(move || operation(&service)).await
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
