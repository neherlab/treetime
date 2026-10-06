use crate::api::extract::{ApiJson, ApiPath};
use crate::api::response::Accepted;
use crate::error::AppError;
use crate::routes::RunPath;
use crate::state::{AppState, LocalSettings};
use aide::axum::ApiRouter;
use aide::axum::routing::{get_with, post_with, put_with};
use app_commands::app_settings::settings::{AnalysisSettings, AppSettings, UiSettings, Workspace, WorkspaceUpdate};
use app_commands::app_settings::workspace::{active_workspace, prepare_workspace};
use app_commands::examples_download::ExamplesDownloadStatus;
use app_commands::runs::record::SaveRunRequest;
use axum::Json;
use axum::extract::State;
use eyre::Report;
use std::sync::Arc;
use treetime_utils::make_internal_report;

pub(crate) fn host_routes() -> ApiRouter<Arc<AppState>> {
  ApiRouter::new().api_route(
    "/api/runs/{id}/save",
    post_with(runs_save, |op| {
      op.id("runsSave").description(
        "Write an output file of a run, or the zip archive of all its outputs, to an absolute path. Only the \
         process that hosts a local app reaches this path; its windows do not.",
      )
    }),
  )
}

pub(crate) fn app_settings_routes() -> ApiRouter<Arc<AppState>> {
  ApiRouter::new()
    .api_route(
      "/api/app-settings",
      get_with(app_settings, |op| {
        op.id("appSettings").description(
          "Settings of a local installation: the folders of the app, the preferences of the user interface, and \
           the settings of the analyses. Only local apps serve this path.",
        )
      }),
    )
    .api_route(
      "/api/app-settings/ui",
      put_with(app_settings_ui, |op| {
        op.id("appSettingsUi")
          .description("Replace the preferences of the user interface. Only local apps serve this path.")
      }),
    )
    .api_route(
      "/api/app-settings/analysis",
      put_with(app_settings_analysis, |op| {
        op.id("appSettingsAnalysis").description(
          "Replace the settings of the analyses; they apply to runs whose configuration does not set them. Only \
           local apps serve this path.",
        )
      }),
    )
    .api_route(
      "/api/examples/download",
      get_with(examples_download, |op| {
        op.id("examplesDownload").description(
          "The download of the example datasets into the examples folder, with the app event that reported it \
           last. Only local apps serve this path.",
        )
      })
      .post_with(examples_download_start, |op| {
        op.id("examplesDownloadStart")
          .description(
            "Start downloading the example datasets into the examples folder; the app events report the progress. \
             Only local apps serve this path.",
          )
          .response_with::<409, AppError, _>(|response| {
            response.description("A download runs already, or the examples folder is not empty")
          })
      }),
    )
    .api_route(
      "/api/workspace",
      get_with(workspace, |op| {
        op.id("workspace")
          .description("The runs folder in use and the default runs folder. Only local apps serve this path.")
      })
      .put_with(workspace_update, |op| {
        op.id("workspaceUpdate").description(
          "Set the runs folder; the folder is created when missing. The back end uses it after it starts again. \
           Only local apps serve this path.",
        )
      }),
    )
}

async fn app_settings(State(state): State<Arc<AppState>>) -> Result<Json<AppSettings>, AppError> {
  blocking(&state, |local| local.store.read()).await.map(Json)
}

async fn app_settings_ui(
  State(state): State<Arc<AppState>>,
  ApiJson(ui): ApiJson<UiSettings>,
) -> Result<Json<UiSettings>, AppError> {
  blocking(&state, move |local| {
    Ok(local.store.update(|settings| settings.ui = ui)?.ui)
  })
  .await
  .map(Json)
}

async fn app_settings_analysis(
  State(state): State<Arc<AppState>>,
  ApiJson(analysis): ApiJson<AnalysisSettings>,
) -> Result<Json<AnalysisSettings>, AppError> {
  blocking(&state, move |local| {
    Ok(local.store.update(|settings| settings.analysis = analysis)?.analysis)
  })
  .await
  .map(Json)
}

async fn workspace(State(state): State<Arc<AppState>>) -> Result<Json<Workspace>, AppError> {
  let local = local_settings(&state)?;
  Ok(Json(active_workspace(&local.paths, local.runs_error.as_deref())))
}

async fn workspace_update(
  State(state): State<Arc<AppState>>,
  ApiJson(WorkspaceUpdate { path }): ApiJson<WorkspaceUpdate>,
) -> Result<Json<AppSettings>, AppError> {
  blocking(&state, move |local| {
    let runs = path
      .as_deref()
      .map(|path| prepare_workspace(&local.paths, path))
      .transpose()?;
    local.store.update(|settings| settings.paths.runs = runs)
  })
  .await
  .map(Json)
}

async fn examples_download(State(state): State<Arc<AppState>>) -> Result<Json<ExamplesDownloadStatus>, AppError> {
  Ok(Json(local_settings(&state)?.examples.status()))
}

async fn examples_download_start(
  State(state): State<Arc<AppState>>,
) -> Result<Accepted<ExamplesDownloadStatus>, AppError> {
  let runs = Arc::clone(&state.runs);
  blocking(&state, move |local| local.examples.start(&runs))
    .await
    .map(Accepted)
}

async fn runs_save(
  State(state): State<Arc<AppState>>,
  ApiPath(RunPath { id }): ApiPath<RunPath>,
  ApiJson(request): ApiJson<SaveRunRequest>,
) -> Result<(), AppError> {
  let runs = Arc::clone(&state.runs);
  blocking(&state, move |_| runs.save(&id, &request)).await
}

async fn blocking<T, F>(state: &AppState, operation: F) -> Result<T, AppError>
where
  T: Send + 'static,
  F: FnOnce(&LocalSettings) -> Result<T, Report> + Send + 'static,
{
  let local = local_settings(state)?;
  Ok(tokio::task::spawn_blocking(move || operation(&local)).await??)
}

fn local_settings(state: &AppState) -> Result<LocalSettings, AppError> {
  Ok(
    state
      .config
      .settings
      .clone()
      .ok_or_else(|| make_internal_report!("the app settings routes were added without app settings"))?,
  )
}
