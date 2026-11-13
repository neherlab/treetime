#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{
    app_with_upload_limit, request, timetree_config, wait_for_status,
  };
  use axum::body::{Body, to_bytes};
  use axum::http::Request;
  use helpers::{Policy, ancestral_package_with_translations, app_with_policy, import_in_chunks, staged_folders};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use tempfile::tempdir;

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_run_package_routes_variant_of_an_imported_run_reads_its_translations_after_the_source_is_deleted_on_desktop()
   {
    check_variant_after_source_deletion(Policy::Desktop).await;
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_run_package_routes_variant_of_an_imported_run_reads_its_translations_after_the_source_is_deleted_in_web2()
   {
    check_variant_after_source_deletion(Policy::Web2).await;
  }

  async fn check_variant_after_source_deletion(policy: Policy) {
    let source = tempdir().unwrap();
    let package = ancestral_package_with_translations(source.path());
    let test = app_with_policy(policy);
    let (imported_status, imported) = import_in_chunks(&test, &package, 4096).await;
    let id = imported["id"].as_str().unwrap().to_owned();
    let (_, settings) = request(&test, "GET", &format!("/api/runs/{id}/settings"), None).await;
    let (created_status, created) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "ancestral", "config": settings, "defer_start": true, "source_run_id": id })),
    )
    .await;
    assert_eq!(200, created_status, "{created}");
    let variant = created["id"].as_str().unwrap().to_owned();
    let (deleted_status, _) = request(&test, "DELETE", &format!("/api/runs/{id}"), None).await;

    let (started_status, started) =
      request(&test, "POST", &format!("/api/runs/{variant}/start"), Some(json!({}))).await;
    let record = wait_for_status(&test, &variant, "ok").await;

    assert_eq!(
      (200, 200, 200, json!("inputs/translations/{cds}.fasta.xz"), 2),
      (
        imported_status,
        deleted_status,
        started_status,
        record["config"]["translations"].clone(),
        record["inputs"]
          .as_array()
          .unwrap()
          .iter()
          .filter(|input| input["setting"] == "translations")
          .count()
      ),
      "{started}"
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_run_package_routes_import_of_a_body_over_the_upload_limit_is_refused_and_leaves_nothing() {
    let test = app_with_upload_limit(1000);
    let (status, body) = import_in_chunks(&test, &[0_u8; 3000], 400).await;
    assert_eq!(
      (
        413,
        json!({
          "code": "upload_too_large",
          "message": "a run package is limited to 1000 bytes, the upload limit, both as uploaded and unpacked; the package is larger",
          "causes": [],
        }),
        0
      ),
      (status, body, staged_folders(&test))
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_run_package_routes_package_body_fails_when_a_packaged_file_disappears() {
    let test = app_with_upload_limit(1 << 30);
    let (_, record) = test.run_with_data("timetree", timetree_config()).await;
    let id = record["id"].as_str().unwrap().to_owned();
    let record = wait_for_status(&test, &id, "ok").await;
    let response = test
      .send(
        Request::get(format!("/api/runs/{id}/package"))
          .body(Body::empty())
          .unwrap(),
      )
      .await;
    let removed = test
      .runs_dir
      .path()
      .canonicalize()
      .unwrap()
      .join(&id)
      .join("out")
      .join(record["output_files"][0]["path"].as_str().unwrap());
    fs::remove_file(&removed).unwrap();

    let body = to_bytes(response.into_body(), usize::MAX).await;

    let Err(error) = body else {
      panic!("expected the body of the run package to fail");
    };
    assert_eq!(format!("When opening '{}'", removed.display()), error.to_string());
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_run_package_routes_import_of_a_package_sent_in_chunks_creates_the_run() {
    let test = app_with_upload_limit(1 << 30);
    let (_, record) = test.run_with_data("timetree", timetree_config()).await;
    let id = record["id"].as_str().unwrap().to_owned();
    wait_for_status(&test, &id, "ok").await;
    let response = test
      .send(
        Request::get(format!("/api/runs/{id}/package"))
          .body(Body::empty())
          .unwrap(),
      )
      .await;
    let package = to_bytes(response.into_body(), usize::MAX).await.unwrap();

    let (status, imported) = import_in_chunks(&test, &package, 1000).await;

    assert_eq!(
      (200, json!("ok"), true, 0),
      (
        status,
        imported["status"].clone(),
        package.len() > 3 * 1000,
        staged_folders(&test)
      )
    );
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::{TestApp, body_json, data_dir};
    use crate::routes::api_router;
    use crate::state::{RouterConfig, confined_service};
    use app_backend::bridge::executor::ThreadExecutor;
    use app_backend::bridge::service::{AppService, Unconfined};
    use app_backend::run_limits::RunLimits;
    use app_backend::runs::manager::RunManager;
    use app_backend::runs::record::{CreateRunRequest, StartRunRequest};
    use app_commands::command::AppCommand;
    use app_commands::json_value::SparseConfig;
    use axum::body::{Body, Bytes};
    use axum::http::Request;
    use serde_json::{Value, json};
    use std::convert::Infallible;
    use std::fs;
    use std::path::Path;
    use std::sync::Arc;
    use tempfile::tempdir;
    use tokio_stream::iter;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Policy {
      Desktop,
      Web2,
    }

    pub(super) fn app_with_policy(policy: Policy) -> TestApp {
      let runs_dir = tempdir().unwrap();
      let service = match policy {
        Policy::Desktop => Arc::new(AppService::new(
          RunManager::open(runs_dir.path()).unwrap(),
          Arc::new(Unconfined),
          RunLimits::Unset,
          Arc::new(ThreadExecutor),
        )),
        Policy::Web2 => confined_service(runs_dir.path(), Arc::new(ThreadExecutor)).unwrap(),
      };
      let (router, _) = api_router(service, RouterConfig::default()).unwrap();
      TestApp { router, runs_dir }
    }

    pub(super) fn ancestral_package_with_translations(root: &Path) -> Vec<u8> {
      let runs = RunManager::open(root).unwrap();
      let data = data_dir().join("rsv/a/20").canonicalize().unwrap();
      let Value::Object(config) = json!({
        "tree": data.join("tree.nwk"),
        "alignment": [data.join("aln.fasta.xz")],
        "translations": data.join("translations/{cds}.fasta.xz"),
        "cds_selection": ["F", "G"],
      }) else {
        unreachable!("the configuration is a mapping");
      };
      let created = runs
        .create(&CreateRunRequest {
          command: AppCommand::Ancestral,
          config: SparseConfig(config),
          defer_start: true,
          source_run_id: None,
          example: None,
        })
        .unwrap();
      runs
        .start(&created.id, &StartRunRequest::default(), Box::new(|_| Ok(())))
        .unwrap()
        .run();
      let record = runs.get(&created.id).unwrap();
      assert_eq!("ok", record.status.to_string(), "{:?}", record.error);
      let mut bytes = vec![];
      runs.package(&created.id).unwrap().stream(&mut bytes).unwrap();
      bytes
    }

    pub(super) async fn import_in_chunks(test: &TestApp, package: &[u8], size: usize) -> (u16, Value) {
      let chunks = package
        .chunks(size)
        .map(|chunk| Ok::<_, Infallible>(Bytes::copy_from_slice(chunk)))
        .collect::<Vec<_>>();
      let response = test
        .send(
          Request::post("/api/runs/import")
            .header("content-type", "application/octet-stream")
            .body(Body::from_stream(iter(chunks)))
            .unwrap(),
        )
        .await;
      body_json(response).await
    }

    pub(super) fn staged_folders(test: &TestApp) -> usize {
      let staging = test.runs_dir.path().join(".staging");
      if staging.exists() {
        fs::read_dir(staging).unwrap().count()
      } else {
        0
      }
    }
  }
}
