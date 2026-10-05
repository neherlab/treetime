#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, request};
  use helpers::{finished_run, local_host};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::fs;
  use std::io::{Cursor, Read};
  use tempfile::tempdir;
  use zip::ZipArchive;

  #[tokio::test]
  async fn test_save_route_writes_an_output_file_to_the_destination() {
    let host = local_host();
    let id = finished_run(&host.app).await;
    let target = tempdir().unwrap();
    let destination = target.path().join("tree.nwk");
    let (status, body) = request(
      &host.app,
      "POST",
      &format!("/api/runs/{id}/save"),
      Some(json!({ "path": "ancestral.nwk", "destination": destination })),
    )
    .await;
    let source = host.app.runs_dir.path().join(&id).join("out").join("ancestral.nwk");
    assert_eq!(
      (200, Value::Null, fs::read(source).unwrap()),
      (status, body, fs::read(&destination).unwrap())
    );
  }

  #[tokio::test]
  async fn test_save_route_writes_the_archive_of_all_outputs_without_a_path() {
    let host = local_host();
    let id = finished_run(&host.app).await;
    let target = tempdir().unwrap();
    let destination = target.path().join("run.zip");
    let (status, _) = request(
      &host.app,
      "POST",
      &format!("/api/runs/{id}/save"),
      Some(json!({ "destination": destination })),
    )
    .await;
    let mut archive = ZipArchive::new(Cursor::new(fs::read(&destination).unwrap())).unwrap();
    let mut tree = String::new();
    archive
      .by_name(&format!("{id}/ancestral.nwk"))
      .unwrap()
      .read_to_string(&mut tree)
      .unwrap();
    let source = fs::read_to_string(host.app.runs_dir.path().join(&id).join("out").join("ancestral.nwk")).unwrap();
    assert_eq!((200, source), (status, tree));
  }

  #[tokio::test]
  async fn test_save_route_answers_an_unknown_run_with_not_found() {
    let host = local_host();
    let target = tempdir().unwrap();
    let (status, error) = request(
      &host.app,
      "POST",
      "/api/runs/missing/save",
      Some(json!({ "destination": target.path().join("run.zip") })),
    )
    .await;
    assert_eq!(
      (404, json!("not_found"), 0),
      (
        status,
        error["code"].clone(),
        fs::read_dir(target.path()).unwrap().count()
      )
    );
  }

  #[tokio::test]
  async fn test_save_route_refuses_a_path_outside_the_run_folder() {
    let host = local_host();
    let id = finished_run(&host.app).await;
    let target = tempdir().unwrap();
    let (status, error) = request(
      &host.app,
      "POST",
      &format!("/api/runs/{id}/save"),
      Some(json!({ "path": "../run.json", "destination": target.path().join("run.json") })),
    )
    .await;
    assert_eq!(
      (
        400,
        json!("file path `../run.json` must name a file inside the run's output folder"),
        0
      ),
      (
        status,
        error["message"].clone(),
        fs::read_dir(target.path()).unwrap().count()
      )
    );
  }

  #[tokio::test]
  async fn test_save_route_refuses_a_relative_destination() {
    let host = local_host();
    let id = finished_run(&host.app).await;
    let (status, error) = request(
      &host.app,
      "POST",
      &format!("/api/runs/{id}/save"),
      Some(json!({ "destination": "run.zip" })),
    )
    .await;
    assert_eq!(
      (
        400,
        json!({
          "code": "invalid_request",
          "message": "the destination must be an absolute path, got 'run.zip'",
          "causes": []
        })
      ),
      (status, error)
    );
  }

  #[tokio::test]
  async fn test_save_route_is_absent_from_a_server_without_local_settings() {
    let server = app();
    let target = tempdir().unwrap();
    let (status, error) = request(
      &server,
      "POST",
      "/api/runs/missing/save",
      Some(json!({ "destination": target.path().join("run.zip") })),
    )
    .await;
    assert_eq!((404, json!("not_found")), (status, error["code"].clone()));
  }

  #[tokio::test]
  async fn test_save_route_is_absent_from_the_renderer_router_of_a_local_app() {
    let renderer = helpers::local_renderer();
    let target = tempdir().unwrap();
    let (status, error) = request(
      &renderer.app,
      "POST",
      "/api/runs/missing/save",
      Some(json!({ "destination": target.path().join("run.zip") })),
    )
    .await;
    assert_eq!((404, json!("not_found")), (status, error["code"].clone()));
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::{TestApp, events_of, request};
    use crate::routes::{LocalRouters, local_api_routers};
    use crate::state::{DEFAULT_MAX_UPLOAD_SIZE, LocalSettings, ServerConfig, server_service};
    use app_commands::app_paths::{AppFolderEnv, AppPaths};
    use app_commands::app_settings::settings::AppPathSettings;
    use app_commands::app_settings::store::AppSettingsStore;
    use serde_json::json;
    use std::sync::Arc;
    use tempfile::{TempDir, tempdir};
    use tokio_util::sync::CancellationToken;

    pub(super) struct LocalApp {
      pub app: TestApp,
      pub _settings: TempDir,
    }

    pub(super) fn local_host() -> LocalApp {
      let (routers, runs_dir, settings) = local_routers();
      LocalApp {
        app: TestApp {
          router: routers.host,
          runs_dir,
        },
        _settings: settings,
      }
    }

    pub(super) fn local_renderer() -> LocalApp {
      let (routers, runs_dir, settings) = local_routers();
      LocalApp {
        app: TestApp {
          router: routers.renderer,
          runs_dir,
        },
        _settings: settings,
      }
    }

    pub(super) async fn finished_run(test: &TestApp) -> String {
      let (_, record) = request(
        test,
        "POST",
        "/api/runs",
        Some(json!({
          "command": "ancestral",
          "config": { "tree": "zika/20/tree.nwk", "alignment": ["zika/20/aln.fasta.xz"] },
        })),
      )
      .await;
      let id = record["id"].as_str().unwrap().to_owned();
      events_of(test, &id, "").await;
      id
    }

    fn local_routers() -> (LocalRouters, TempDir, TempDir) {
      let dir = tempdir().unwrap();
      let runs_dir = tempdir().unwrap();
      let paths = AppPaths::resolve(
        dir.path(),
        &AppFolderEnv::default(),
        &AppPathSettings {
          runs: Some(runs_dir.path().to_path_buf()),
          ..AppPathSettings::default()
        },
      );
      let config = ServerConfig {
        examples_dir: TestApp::examples_dir(),
        runs_dir: paths.runs.path.clone(),
        max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
        shutdown: CancellationToken::new(),
        settings: Some(LocalSettings {
          store: Arc::new(AppSettingsStore::open(dir.path()).unwrap()),
          paths,
          runs_error: None,
        }),
      };
      let routers = local_api_routers(server_service(&config).unwrap(), config).unwrap();
      (routers, runs_dir, dir)
    }
  }
}
