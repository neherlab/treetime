use crate::api::extract::ApiJson;
use crate::error::AppError;
use crate::state::{AppState, LocalSettings};
use aide::axum::ApiRouter;
use aide::axum::routing::{get_with, put_with};
use app_commands::app_settings::settings::{AppSettings, UiSettings, Workspace, WorkspaceUpdate};
use app_commands::app_settings::workspace::prepare_workspace;
use axum::Json;
use axum::extract::State;
use eyre::Report;
use std::sync::Arc;
use treetime_utils::make_internal_report;

pub(crate) fn app_settings_routes() -> ApiRouter<Arc<AppState>> {
  ApiRouter::new()
    .api_route(
      "/api/app-settings",
      get_with(app_settings, |op| {
        op.id("appSettings").description(
          "Settings of a local installation: the runs folder and the preferences of the user interface. Only \
           local apps serve this path.",
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

async fn workspace(State(state): State<Arc<AppState>>) -> Result<Json<Workspace>, AppError> {
  let local = local_settings(&state)?;
  Ok(Json(Workspace {
    path: state.config.runs_dir.clone(),
    default_path: local.default_workspace,
  }))
}

async fn workspace_update(
  State(state): State<Arc<AppState>>,
  ApiJson(WorkspaceUpdate { path }): ApiJson<WorkspaceUpdate>,
) -> Result<Json<AppSettings>, AppError> {
  blocking(&state, move |local| {
    let workspace = path.as_deref().map(prepare_workspace).transpose()?;
    local.store.update(|settings| settings.workspace = workspace)
  })
  .await
  .map(Json)
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
