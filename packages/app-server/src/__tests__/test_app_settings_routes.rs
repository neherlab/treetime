#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, request};
  use app_commands::app_paths::AppFolderEnv;
  use helpers::{local_app, local_app_with};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use std::fs;
  use tempfile::tempdir;

  #[tokio::test]
  async fn test_app_settings_routes_read_defaults_without_a_settings_file() {
    let local = local_app();
    let actual = request(&local.app, "GET", "/api/app-settings", None).await;
    assert_eq!((200, json!({})), actual);
  }

  #[tokio::test]
  async fn test_app_settings_routes_write_the_ui_preferences_to_the_settings_file() {
    let local = local_app();
    let body = json!({ "theme": "dark", "sidebar_width": 420 });
    let (status, ui) = request(&local.app, "PUT", "/api/app-settings/ui", Some(body.clone())).await;
    let expected_file = indoc! {r#"
      ui:
        theme: "dark"
        sidebar_width: 420
    "#};
    assert_eq!(
      (200, body, expected_file.to_owned()),
      (
        status,
        ui,
        fs::read_to_string(local.dir.path().join("settings.yaml")).unwrap()
      )
    );
  }

  #[tokio::test]
  async fn test_app_settings_routes_refuse_an_unknown_ui_preference() {
    let local = local_app();
    let (status, error) = request(
      &local.app,
      "PUT",
      "/api/app-settings/ui",
      Some(json!({ "color": "red" })),
    )
    .await;
    assert_eq!((400, json!("invalid_request")), (status, error["code"].clone()));
  }

  #[tokio::test]
  async fn test_app_settings_routes_report_the_active_and_default_workspace() {
    let local = local_app();
    let runs = local.app.runs_dir.path().to_path_buf();
    let actual = request(&local.app, "GET", "/api/workspace", None).await;
    let expected = json!({ "path": runs, "default_path": local.dir.path().join("runs"), "fixed_by": null });
    assert_eq!((200, expected), actual);
  }

  #[tokio::test]
  async fn test_app_settings_routes_store_a_new_workspace_and_create_it() {
    let local = local_app();
    let folder = local.dir.path().join("elsewhere").join("runs");
    let (status, settings) = request(&local.app, "PUT", "/api/workspace", Some(json!({ "path": folder }))).await;
    assert_eq!(
      (200, json!({ "paths": { "runs": folder } }), true),
      (status, settings, folder.is_dir())
    );
  }

  #[tokio::test]
  async fn test_app_settings_routes_reset_the_workspace_to_the_default() {
    let local = local_app();
    let folder = local.dir.path().join("elsewhere");
    request(&local.app, "PUT", "/api/workspace", Some(json!({ "path": folder }))).await;
    let actual = request(&local.app, "PUT", "/api/workspace", Some(json!({ "path": null }))).await;
    assert_eq!((200, json!({})), actual);
  }

  #[tokio::test]
  async fn test_app_settings_routes_refuse_a_relative_workspace() {
    let local = local_app();
    let (status, error) = request(&local.app, "PUT", "/api/workspace", Some(json!({ "path": "runs" }))).await;
    assert_eq!(
      (
        400,
        json!({
          "code": "invalid_request",
          "message": "the runs folder must be an absolute path, got 'runs'",
          "causes": [],
        })
      ),
      (status, error)
    );
  }

  #[tokio::test]
  async fn test_app_settings_routes_refuse_a_workspace_fixed_by_the_environment() {
    let scratch = tempdir().unwrap();
    let local = local_app_with(&AppFolderEnv {
      runs: Some(scratch.path().to_path_buf()),
      ..AppFolderEnv::default()
    });
    let (status, error) = request(
      &local.app,
      "PUT",
      "/api/workspace",
      Some(json!({ "path": "/data/runs" })),
    )
    .await;
    assert_eq!(
      (
        400,
        json!(
          "the environment variable TREETIME_RUNS_DIR sets the runs folder; unset it to choose the folder in the app"
        )
      ),
      (status, error["message"].clone())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::settings(          "GET", "/api/app-settings")]
  #[case::ui(                "PUT", "/api/app-settings/ui")]
  #[case::workspace(         "GET", "/api/workspace")]
  #[case::workspace_update(  "PUT", "/api/workspace")]
  #[trace]
  #[tokio::test]
  async fn test_app_settings_routes_are_absent_from_a_server_without_local_settings(
    #[case] method: &str,
    #[case] uri: &str,
  ) {
    let server = app();
    let (status, error) = request(&server, method, uri, Some(json!({}))).await;
    assert_eq!((404, json!("not_found")), (status, error["code"].clone()));
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::TestApp;
    use crate::create_router;
    use crate::state::{DEFAULT_MAX_UPLOAD_SIZE, LocalSettings, ServerConfig};
    use crate::web::WebOptions;
    use app_commands::app_paths::{AppFolderEnv, AppPaths};
    use app_commands::app_settings::settings::AppPathSettings;
    use app_commands::app_settings::store::AppSettingsStore;
    use std::sync::Arc;
    use tempfile::{TempDir, tempdir};
    use tokio_util::sync::CancellationToken;

    pub(super) struct LocalApp {
      pub app: TestApp,
      pub dir: TempDir,
    }

    pub(super) fn local_app() -> LocalApp {
      local_app_with(&AppFolderEnv::default())
    }

    pub(super) fn local_app_with(env: &AppFolderEnv) -> LocalApp {
      let dir = tempdir().unwrap();
      let runs_dir = tempdir().unwrap();
      let paths = AppPaths::resolve(
        dir.path(),
        env,
        &AppPathSettings {
          runs: Some(runs_dir.path().to_path_buf()),
          ..AppPathSettings::default()
        },
      );
      let router = create_router(
        ServerConfig {
          data_dir: TestApp::data_dir(),
          runs_dir: paths.runs.path.clone(),
          max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
          shutdown: CancellationToken::new(),
          settings: Some(LocalSettings {
            store: Arc::new(AppSettingsStore::open(dir.path()).unwrap()),
            paths,
          }),
        },
        &WebOptions::default(),
      )
      .unwrap();
      LocalApp {
        app: TestApp { router, runs_dir },
        dir,
      }
    }
  }
}
