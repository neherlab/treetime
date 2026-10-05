#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, open_app_events, request, take_events};
  use helpers::{local_app_with_examples, serve_examples};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::fs;

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_examples_routes_download_into_the_examples_folder_and_report_it_on_the_app_events() {
    let server = serve_examples();
    let local = local_app_with_examples(&server.url("examples.zip"));
    let mut events = open_app_events(&local.app, "", None).await;

    let (status, started) = request(&local.app, "POST", "/api/examples/download", None).await;
    server.release();
    let mut download_events = vec![];
    while download_events
      .last()
      .is_none_or(|data: &Value| data["state"] == "running")
    {
      let event = take_events(&mut events, 1).await.remove(0);
      download_events.push(event.data["download"].clone());
      assert_eq!(json!("examples-download"), event.data["kind"]);
    }
    let last = download_events.last().unwrap().clone();
    let (_, state) = request(&local.app, "GET", "/api/examples/download", None).await;
    let (_, catalog) = request(&local.app, "GET", "/api/datasets", None).await;

    assert_eq!(
      (
        202,
        json!({ "state": "running", "received": 0 }),
        json!("done"),
        json!("done"),
        json!(["zika/20"])
      ),
      (
        status,
        started["download"].clone(),
        last["state"].clone(),
        state["download"]["state"].clone(),
        json!(
          catalog["datasets"]
            .as_array()
            .unwrap()
            .iter()
            .map(|dataset| dataset["name"].clone())
            .collect::<Vec<_>>()
        )
      )
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_examples_routes_state_names_the_event_that_reported_it_last() {
    let server = serve_examples();
    let local = local_app_with_examples(&server.url("examples.zip"));
    let mut events = open_app_events(&local.app, "", None).await;
    request(&local.app, "POST", "/api/examples/download", None).await;
    server.release();
    let mut last = take_events(&mut events, 1).await.remove(0);
    while last.data["download"]["state"] == "running" {
      last = take_events(&mut events, 1).await.remove(0);
    }
    let (status, state) = request(&local.app, "GET", "/api/examples/download", None).await;
    assert_eq!(
      (200, json!(last.data["seq"]), last.data["download"].clone()),
      (status, state["seq"].clone(), state["download"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_examples_routes_refuse_a_second_download_while_one_runs() {
    let server = serve_examples();
    let local = local_app_with_examples(&server.url("examples.zip"));
    let (first, _) = request(&local.app, "POST", "/api/examples/download", None).await;
    let (second, error) = request(&local.app, "POST", "/api/examples/download", None).await;
    server.release();
    assert_eq!(
      (
        202,
        409,
        json!({ "code": "conflict", "message": "the example datasets are being downloaded already", "causes": [] })
      ),
      (first, second, error)
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_examples_routes_refuse_a_download_into_a_folder_that_is_not_empty() {
    let server = serve_examples();
    let local = local_app_with_examples(&server.url("examples.zip"));
    let examples = local.dir.path().join("examples");
    fs::create_dir_all(&examples).unwrap();
    fs::write(examples.join("mine.nwk"), "(A:1);\n").unwrap();
    let (status, error) = request(&local.app, "POST", "/api/examples/download", None).await;
    let (_, state) = request(&local.app, "GET", "/api/examples/download", None).await;
    assert_eq!(
      (
        409,
        json!(format!(
          "the folder '{}' is not empty; the example datasets go into an empty or missing folder",
          examples.display()
        )),
        json!({ "download": { "state": "idle" } })
      ),
      (status, error["message"].clone(), state)
    );
  }

  #[tokio::test]
  async fn test_examples_routes_are_absent_from_a_server_without_local_settings() {
    let server = app();
    let (status, error) = request(&server, "POST", "/api/examples/download", None).await;
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
    use app_commands::examples_download::ExampleDownloads;
    use axum::Router;
    use axum::routing::get;
    use std::fs;
    use std::io::{Cursor, Write};
    use std::net::SocketAddr;
    use std::sync::{Arc, mpsc};
    use std::thread;
    use tempfile::{TempDir, tempdir};
    use tokio::sync::Notify;
    use tokio_util::sync::CancellationToken;
    use zip::ZipWriter;
    use zip::write::SimpleFileOptions;

    pub(super) struct LocalApp {
      pub app: TestApp,
      pub dir: TempDir,
    }

    pub(super) struct ExamplesServer {
      address: SocketAddr,
      gate: Arc<Notify>,
    }

    impl ExamplesServer {
      pub(super) fn url(&self, name: &str) -> String {
        format!("http://{}/{name}", self.address)
      }

      pub(super) fn release(&self) {
        self.gate.notify_one();
      }
    }

    pub(super) fn serve_examples() -> ExamplesServer {
      let gate = Arc::new(Notify::new());
      let waiting = Arc::clone(&gate);
      let (send, receive) = mpsc::sync_channel(1);
      thread::spawn(move || {
        let runtime = tokio::runtime::Builder::new_current_thread()
          .enable_all()
          .build()
          .unwrap();
        runtime.block_on(async move {
          let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
          send.send(listener.local_addr().unwrap()).unwrap();
          let router = Router::new().route(
            "/examples.zip",
            get(move || {
              let waiting = Arc::clone(&waiting);
              async move {
                waiting.notified().await;
                archive()
              }
            }),
          );
          axum::serve(listener, router).await.unwrap();
        });
      });
      ExamplesServer {
        address: receive.recv().unwrap(),
        gate,
      }
    }

    pub(super) fn local_app_with_examples(url: &str) -> LocalApp {
      let dir = tempdir().unwrap();
      let runs_dir = tempdir().unwrap();
      let examples = dir.path().join("examples");
      fs::create_dir_all(&examples).unwrap();
      let paths = AppPaths::resolve(
        dir.path(),
        &AppFolderEnv::default(),
        &AppPathSettings {
          runs: Some(runs_dir.path().to_path_buf()),
          ..AppPathSettings::default()
        },
      );
      let router = create_router(
        ServerConfig {
          examples_dir: examples.clone(),
          runs_dir: paths.runs.path.clone(),
          max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
          shutdown: CancellationToken::new(),
          settings: Some(LocalSettings {
            store: Arc::new(AppSettingsStore::open(dir.path()).unwrap()),
            paths,
            runs_error: None,
            examples: Arc::new(ExampleDownloads::new(url, &examples)),
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

    fn archive() -> Vec<u8> {
      let mut writer = ZipWriter::new(Cursor::new(Vec::new()));
      writer
        .start_file("zika/20/tree.nwk", SimpleFileOptions::default())
        .unwrap();
      writer.write_all(b"(A:1,B:1);\n").unwrap();
      writer.finish().unwrap().into_inner()
    }
  }
}
