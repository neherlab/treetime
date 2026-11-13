#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, request};
  use app_backend::app_paths::AppFolderEnv;
  use helpers::{local_app, local_app_with};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
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
  async fn test_app_settings_routes_write_the_analysis_settings_to_the_settings_file() {
    let local = local_app();
    let body = json!({ "max_grid_points": 250000 });
    let (status, analysis) = request(&local.app, "PUT", "/api/app-settings/analysis", Some(body.clone())).await;
    let (_, settings) = request(&local.app, "GET", "/api/app-settings", None).await;
    let expected_file = indoc! {r#"
      analysis:
        max_grid_points: 250000
    "#};
    assert_eq!(
      (200, body.clone(), json!({ "analysis": body }), expected_file.to_owned()),
      (
        status,
        analysis,
        settings,
        fs::read_to_string(local.dir.path().join("settings.yaml")).unwrap()
      )
    );
  }

  #[tokio::test]
  async fn test_app_settings_routes_clear_the_analysis_settings() {
    let local = local_app();
    request(
      &local.app,
      "PUT",
      "/api/app-settings/analysis",
      Some(json!({ "max_grid_points": 250000 })),
    )
    .await;
    let actual = request(&local.app, "PUT", "/api/app-settings/analysis", Some(json!({}))).await;
    let (_, settings) = request(&local.app, "GET", "/api/app-settings", None).await;
    assert_eq!(((200, json!({})), json!({})), (actual, settings));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_minimum(json!({ "max_grid_points": 999 }))]
  #[case::not_a_number( json!({ "max_grid_points": "many" }))]
  #[case::unknown_key(  json!({ "grid": 1000 }))]
  #[trace]
  #[tokio::test]
  async fn test_app_settings_routes_refuse_invalid_analysis_settings(#[case] body: Value) {
    let local = local_app();
    let (status, error) = request(&local.app, "PUT", "/api/app-settings/analysis", Some(body)).await;
    assert_eq!(
      (400, json!("invalid_request"), false),
      (status, error["code"].clone(), local.dir.path().join("settings.yaml").exists())
    );
  }

  #[tokio::test]
  async fn test_app_settings_routes_report_the_active_and_default_workspace() {
    let local = local_app();
    let runs = local.app.runs_dir.path().to_path_buf();
    let actual = request(&local.app, "GET", "/api/workspace", None).await;
    let expected = json!({ "path": runs, "default_path": local.dir.path().join("runs") });
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
    let actual = request(&local.app, "PUT", "/api/workspace", Some(json!({}))).await;
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
  #[case::analysis(          "PUT", "/api/app-settings/analysis")]
  #[case::workspace(         "GET", "/api/workspace")]
  #[case::workspace_update(  "PUT", "/api/workspace")]
  #[trace]
  #[tokio::test]
  async fn test_app_settings_routes_are_absent_from_the_public_router(
    #[case] method: &str,
    #[case] uri: &str,
  ) {
    let public = app();
    let (status, error) = request(&public, method, uri, Some(json!({}))).await;
    assert_eq!((404, json!("not_found")), (status, error["code"].clone()));
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::{TestApp, local_routers};
    use app_backend::app_paths::AppFolderEnv;
    use app_backend::examples_download::ExampleDownloads;
    use tempfile::{TempDir, tempdir};

    const EXAMPLES_URL: &str = "http://127.0.0.1:9/examples.zip";

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
      let router = local_routers(
        dir.path(),
        runs_dir.path(),
        env,
        ExampleDownloads::new(EXAMPLES_URL, dir.path().join("examples")),
      )
      .renderer;
      LocalApp {
        app: TestApp { router, runs_dir },
        dir,
      }
    }
  }
}
