#[cfg(test)]
mod tests {
  use crate::routes::api_doc;
  use helpers::{
    TestApp, app, app_with_upload_limit, events_of, next_event, read_events, request, timetree_config, wait_for_status,
  };
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_run_streams_events_and_records_outputs() {
    let test = app();
    let (status, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "timetree", "config": timetree_config(), "title": "zika" })),
    )
    .await;
    assert_eq!(
      (200, json!("running"), json!("zika")),
      (status, record["status"].clone(), record["title"].clone())
    );
    let id = record["id"].as_str().unwrap().to_owned();

    let events = events_of(&test, &id, "").await;
    let types: Vec<&str> = events.iter().map(|event| event.name.as_str()).collect();
    assert_eq!(Some(&"started"), types.first());
    assert_eq!(Some(&"terminal"), types.last());
    assert_eq!(1, types.iter().filter(|name| **name == "terminal").count());
    assert!(types.contains(&"progress") && types.contains(&"iteration"), "{types:?}");
    assert!(
      events
        .iter()
        .enumerate()
        .all(|(index, event)| event.data["seq"] == json!(index)),
      "events are numbered from 0 without gaps"
    );

    let (_, record) = request(&test, "GET", &format!("/api/runs/{id}"), None).await;
    let run_dir = test.runs_dir.path().canonicalize().unwrap().join(&id);
    assert_eq!(
      (
        json!("ok"),
        json!(run_dir.join("out")),
        json!([
          "Nwk",
          "Nexus",
          "Auspice",
          "AugurNodeData",
          "Gtr",
          "ReconstructedNucFasta",
          "ClockModel",
          "CoalescentTsv",
          "Tracelog"
        ])
      ),
      (
        record["status"].clone(),
        record["config"]["output_all"].clone(),
        record["config"]["output_selection"].clone()
      )
    );
    assert!(record["headline"]["clock_rate"].as_f64().unwrap() > 0.0);
    assert_eq!(64, record["config_hash"].as_str().unwrap().len());

    let (_, files) = request(&test, "GET", &format!("/api/runs/{id}/files"), None).await;
    let kinds: Vec<(String, Value)> = files
      .as_array()
      .unwrap()
      .iter()
      .map(|file| (file["path"].as_str().unwrap().to_owned(), file["kind"].clone()))
      .filter(|(path, _)| path.ends_with(".auspice.json") || path.ends_with(".tracelog.csv"))
      .collect();
    assert_eq!(
      vec![
        ("timetree.auspice.json".to_owned(), json!("auspice")),
        ("timetree.tracelog.csv".to_owned(), json!("tracelog")),
      ],
      kinds
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_events_resume_from_an_offset() {
    let test = app();
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": { "tree": "zika/20/tree.nwk", "metadata": "zika/20/metadata.tsv" } })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let all = events_of(&test, &id, "").await;
    let from_three = events_of(&test, &id, "?from=3").await;
    let seqs = |events: &[helpers::SseEvent]| events.iter().map(|event| event.data["seq"].clone()).collect::<Vec<_>>();
    assert_eq!(seqs(&all[3..]), seqs(&from_three));

    let response = test
      .send(
        axum::http::Request::get(format!("/api/runs/{id}/events"))
          .header("last-event-id", "4")
          .body(axum::body::Body::empty())
          .unwrap(),
      )
      .await;
    assert_eq!(seqs(&all[5..]), seqs(&read_events(response).await));
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_rejected_configuration_ends_with_an_error_run() {
    let test = app();
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": { "tree": "zika/20/tree.nwk", "bogus": 1 } })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let events = events_of(&test, &id, "").await;
    let terminal = &events.last().unwrap().data;
    assert_eq!(
      (json!("error"), json!("invalid configuration: unknown field `bogus`")),
      (terminal["data"]["status"].clone(), terminal["data"]["message"].clone())
    );
    let record = wait_for_status(&test, &id, "error").await;
    assert_eq!(
      json!("invalid configuration: unknown field `bogus`"),
      record["error"]["message"]
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_input_outside_the_data_dir_ends_with_an_error_run() {
    let test = app();
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": { "tree": "../Cargo.toml", "metadata": "zika/20/metadata.tsv" } })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let events = events_of(&test, &id, "").await;
    assert_eq!(
      json!("input `../Cargo.toml` of setting `tree` is outside the directories the server reads inputs from"),
      events.last().unwrap().data["data"]["message"]
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_cancel_one_of_two_concurrent_runs() {
    let test = app();
    let create = || {
      request(
        &test,
        "POST",
        "/api/runs",
        Some(json!({ "command": "timetree", "config": timetree_config() })),
      )
    };
    let (_, first) = create().await;
    let (_, second) = create().await;
    let first = first["id"].as_str().unwrap().to_owned();
    let second = second["id"].as_str().unwrap().to_owned();

    let mut stream = test
      .send(
        axum::http::Request::get(format!("/api/runs/{first}/events"))
          .body(axum::body::Body::empty())
          .unwrap(),
      )
      .await
      .into_body()
      .into_data_stream();
    let mut buffer = String::new();
    let started = next_event(&mut stream, &mut buffer).await.unwrap();
    assert_eq!("started", started.name);

    let (status, response) = request(&test, "POST", &format!("/api/runs/{first}/cancel"), None).await;
    assert_eq!((200, json!(true)), (status, response["cancelled"].clone()));

    let first_record = wait_for_status(&test, &first, "cancelled").await;
    let second_record = wait_for_status(&test, &second, "ok").await;
    assert_eq!(
      (json!("cancelled"), json!("ok")),
      (first_record["status"].clone(), second_record["status"].clone())
    );
    let second_events = events_of(&test, &second, "").await;
    assert!(second_events.iter().any(|event| event.name == "progress"));

    let (status, response) = request(&test, "POST", &format!("/api/runs/{first}/cancel"), None).await;
    assert_eq!(
      (200, json!(false)),
      (status, response["cancelled"].clone()),
      "a finished run cannot be cancelled"
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_list_rename_pin_delete_restore_and_purge() {
    let test = app();
    let mut ids = vec![];
    for title in ["older", "newer"] {
      let (_, record) = request(
        &test,
        "POST",
        "/api/runs",
        Some(json!({ "command": "clock", "config": {}, "title": title, "defer_start": true })),
      )
      .await;
      ids.push(record["id"].as_str().unwrap().to_owned());
    }
    let (_, list) = request(&test, "GET", "/api/runs", None).await;
    let titles: Vec<&str> = list["runs"]
      .as_array()
      .unwrap()
      .iter()
      .map(|run| run["title"].as_str().unwrap())
      .collect();
    assert_eq!(
      (vec!["newer", "older"], json!(0)),
      (titles, list["active_runs"].clone())
    );

    let (status, summary) = request(
      &test,
      "PATCH",
      &format!("/api/runs/{}", ids[0]),
      Some(json!({ "title": "renamed", "pinned": true })),
    )
    .await;
    assert_eq!(
      (200, json!("renamed"), json!(true)),
      (status, summary["title"].clone(), summary["pinned"].clone())
    );

    let (status, _) = request(&test, "DELETE", &format!("/api/runs/{}", ids[0]), None).await;
    assert_eq!(204, status);
    let (status, _) = request(&test, "GET", &format!("/api/runs/{}", ids[0]), None).await;
    assert_eq!(404, status);
    let (status, summary) = request(&test, "POST", &format!("/api/runs/{}/restore", ids[0]), None).await;
    assert_eq!((200, json!("renamed")), (status, summary["title"].clone()));

    request(&test, "DELETE", &format!("/api/runs/{}", ids[1]), None).await;
    let (status, _) = request(&test, "POST", &format!("/api/runs/{}/purge", ids[1]), None).await;
    assert_eq!(204, status);
    let (status, _) = request(&test, "POST", &format!("/api/runs/{}/restore", ids[1]), None).await;
    assert_eq!(404, status, "a purged run cannot be restored");
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_file_paths_that_leave_the_output_folder_are_rejected() {
    let test = app();
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": { "tree": "zika/20/tree.nwk", "metadata": "zika/20/metadata.tsv" } })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    wait_for_status(&test, &id, "ok").await;

    let file = |path: &str| format!("/api/runs/{id}/file?path={path}");
    let (status, body) = request(&test, "GET", &file("..%2Frun.json"), None).await;
    assert_eq!(
      (
        400,
        json!("file path `../run.json` must name a file inside the run's output folder")
      ),
      (status, body["message"].clone())
    );
    let (status, _) = request(&test, "GET", &file("%2Fetc%2Fpasswd"), None).await;
    assert_eq!(400, status);
    let (status, _) = request(&test, "GET", &file("missing.json"), None).await;
    assert_eq!(404, status);
    let (status, model) = request(&test, "GET", &file("clock.clock-model.json"), None).await;
    assert_eq!((200, true), (status, model["clock_rate"].as_f64().unwrap() > 0.0));

    let response = test
      .send(
        axum::http::Request::get(format!("/api/runs/{id}/archive"))
          .body(axum::body::Body::empty())
          .unwrap(),
      )
      .await;
    let content_type = response.headers()["content-type"].clone();
    let bytes = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
    assert_eq!(
      ("application/zip", &b"PK"[..]),
      (content_type.to_str().unwrap(), &bytes[..2])
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_uploaded_inputs_feed_the_run() {
    let test = app();
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": {}, "defer_start": true })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let data = TestApp::data_dir();
    let tree = test
      .upload(&id, "tree.nwk", std::fs::read(data.join("zika/20/tree.nwk")).unwrap())
      .await;
    let metadata = test
      .upload(
        &id,
        "metadata.tsv",
        std::fs::read(data.join("zika/20/metadata.tsv")).unwrap(),
      )
      .await;
    assert_eq!(200, tree.0);

    let (status, record) = request(
      &test,
      "POST",
      &format!("/api/runs/{id}/start"),
      Some(json!({ "config": { "tree": tree.1["path"], "metadata": metadata.1["path"] } })),
    )
    .await;
    assert_eq!((200, json!("running")), (status, record["status"].clone()));
    let record = wait_for_status(&test, &id, "ok").await;
    assert_eq!(
      json!(
        test
          .runs_dir
          .path()
          .canonicalize()
          .unwrap()
          .join(&id)
          .join("inputs/tree.nwk")
      ),
      record["inputs"][0]["path"]
    );

    let late = test.upload(&id, "late.nwk", b"(A,B);".to_vec()).await;
    assert_eq!(409, late.0, "inputs cannot change after the run starts");
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_upload_over_the_limit_names_the_limit() {
    let test = app_with_upload_limit(1000);
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": {}, "defer_start": true })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let small = test.upload(&id, "a.nwk", vec![b'x'; 600]).await;
    let too_much = test.upload(&id, "b.nwk", vec![b'x'; 600]).await;
    let replace = test.upload(&id, "a.nwk", vec![b'x'; 900]).await;
    assert_eq!(
      (
        200,
        413,
        json!("the inputs of a run are limited to 1000 bytes in total; `b.nwk` does not fit"),
        200
      ),
      (small.0, too_much.0, too_much.1["message"].clone(), replace.0)
    );
    let (status, body) = test.upload(&id, "..", vec![b'x'; 1]).await;
    assert_eq!(400, status, "{body}");
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_check_inputs_reports_facts_of_zika_86() {
    let test = app();
    let (status, facts) = request(
      &test,
      "POST",
      "/api/check-inputs",
      Some(json!({
        "tree": "zika/86/tree.nwk",
        "metadata": "zika/86/metadata.tsv",
        "alignment": ["zika/86/aln.fasta.xz"],
      })),
    )
    .await;
    assert_eq!(
      (
        200,
        json!(86),
        json!(15),
        json!(86),
        json!(54),
        json!("name"),
        json!("date"),
        json!([]),
        json!([]),
        json!(86)
      ),
      (
        status,
        facts["tree"]["tips"].clone(),
        facts["tree"]["polytomies"].clone(),
        facts["metadata"]["dates"]["exact_days"].clone(),
        facts["metadata"]["dates"]["on_day_1_or_15"].clone(),
        facts["metadata"]["id_column"].clone(),
        facts["metadata"]["date_column"].clone(),
        facts["tips_without_metadata"].clone(),
        facts["tips_without_sequence"].clone(),
        facts["alignment"]["sequences"].clone(),
      )
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_check_inputs_rejects_paths_outside_the_data_dir() {
    let test = app();
    let (status, body) = request(
      &test,
      "POST",
      "/api/check-inputs",
      Some(json!({ "tree": "../Cargo.toml" })),
    )
    .await;
    assert_eq!(
      (
        400,
        json!("input `../Cargo.toml` of setting `tree` is outside the directories the server reads inputs from")
      ),
      (status, body["message"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_check_config_reports_cli_errors() {
    let test = app();
    let (_, response) = request(
      &test,
      "POST",
      "/api/check-config",
      Some(json!({ "command": "prune", "text": "tree: t.nwk\nprune_short: fast\n" })),
    )
    .await;
    assert_eq!(
      (
        json!("invalid"),
        json!("invalid configuration: \"fast\" is not of types \"null\", \"number\"")
      ),
      (response["status"].clone(), response["message"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_datasets_lists_examples() {
    let test = app();
    let (_, catalog) = request(&test, "GET", "/api/datasets", None).await;
    let example = catalog["examples"]
      .as_array()
      .unwrap()
      .iter()
      .find(|example| example["path"] == json!("ebola/20/ancestral-parsimony.yaml"))
      .unwrap();
    assert_eq!(json!("ancestral"), example["command"]);
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_invalid_run_id_is_a_bad_request() {
    let test = app();
    let (status, _) = request(&test, "GET", "/api/runs/a.b", None).await;
    assert_eq!(400, status);
  }

  #[test]
  fn test_routes_openapi_describes_run_bodies_and_bridge_types() {
    let doc = api_doc().unwrap();
    assert_eq!(
      (
        json!("request"),
        json!("#/components/schemas/CreateRunRequest"),
        json!("stream"),
        json!("#/components/schemas/RunEvent"),
        json!("upload"),
        json!("query"),
      ),
      (
        doc["paths"]["/api/runs"]["post"]["x-bridge-type"].clone(),
        doc["paths"]["/api/runs"]["post"]["requestBody"]["content"]["application/json"]["schema"]["$ref"].clone(),
        doc["paths"]["/api/runs/{id}/events"]["get"]["x-bridge-type"].clone(),
        doc["paths"]["/api/runs/{id}/events"]["get"]["responses"]["200"]["content"]["text/event-stream"]["schema"]
          ["$ref"]
          .clone(),
        doc["paths"]["/api/runs/{id}/inputs/{name}"]["put"]["x-bridge-type"].clone(),
        doc["paths"]["/api/version"]["get"]["x-bridge-type"].clone(),
      )
    );
    let config = &doc["components"]["schemas"]["ClockConfig"];
    assert_eq!(
      (json!(false), json!("#/components/schemas/BranchSplitArgs")),
      (
        config["additionalProperties"].clone(),
        config["properties"]["branch_split"]["$ref"].clone()
      )
    );
  }

  mod helpers {
    use crate::create_router;
    use crate::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig};
    use axum::Router;
    use axum::body::{Body, BodyDataStream};
    use axum::http::Request;
    use axum::response::Response;
    use serde_json::{Value, json};
    use std::path::PathBuf;
    use std::str;
    use std::time::Duration;
    use tempfile::{TempDir, tempdir};
    use tokio_stream::StreamExt;
    use tower::ServiceExt;

    pub(super) struct SseEvent {
      pub name: String,
      pub data: Value,
    }

    pub(super) struct TestApp {
      pub router: Router,
      pub runs_dir: TempDir,
    }

    impl TestApp {
      pub(super) fn data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data")
      }

      pub(super) async fn send(&self, request: Request<Body>) -> Response {
        self.router.clone().oneshot(request).await.unwrap()
      }

      pub(super) async fn upload(&self, id: &str, name: &str, bytes: Vec<u8>) -> (u16, Value) {
        let response = self
          .send(
            Request::put(format!("/api/runs/{id}/inputs/{name}"))
              .body(Body::from(bytes))
              .unwrap(),
          )
          .await;
        body_json(response).await
      }
    }

    pub(super) fn app() -> TestApp {
      app_with_upload_limit(DEFAULT_MAX_UPLOAD_SIZE)
    }

    pub(super) fn app_with_upload_limit(max_upload_size: usize) -> TestApp {
      let runs_dir = tempdir().unwrap();
      let router = create_router(
        ServerConfig {
          data_dir: TestApp::data_dir(),
          runs_dir: runs_dir.path().to_path_buf(),
          max_upload_size,
        },
        None,
      )
      .unwrap();
      TestApp { router, runs_dir }
    }

    pub(super) fn timetree_config() -> Value {
      json!({
        "tree": "zika/20/tree.nwk",
        "metadata": "zika/20/metadata.tsv",
        "alignment": ["zika/20/aln.fasta.xz"],
        "max_iter": 2,
        "seed": 7,
      })
    }

    pub(super) async fn request(test: &TestApp, method: &str, uri: &str, body: Option<Value>) -> (u16, Value) {
      let builder = Request::builder()
        .method(method)
        .uri(uri)
        .header("content-type", "application/json");
      let request = match body {
        Some(body) => builder.body(Body::from(body.to_string())),
        None => builder.body(Body::empty()),
      }
      .unwrap();
      body_json(test.send(request).await).await
    }

    pub(super) async fn body_json(response: Response) -> (u16, Value) {
      let status = response.status().as_u16();
      let bytes = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
      let value = if bytes.is_empty() {
        Value::Null
      } else {
        serde_json::from_slice(&bytes).unwrap_or_else(|_| json!(String::from_utf8_lossy(&bytes)))
      };
      (status, value)
    }

    pub(super) async fn events_of(test: &TestApp, id: &str, query: &str) -> Vec<SseEvent> {
      let response = test
        .send(
          Request::get(format!("/api/runs/{id}/events{query}"))
            .body(Body::empty())
            .unwrap(),
        )
        .await;
      read_events(response).await
    }

    pub(super) async fn wait_for_status(test: &TestApp, id: &str, status: &str) -> Value {
      for _ in 0..600 {
        let (_, record) = request(test, "GET", &format!("/api/runs/{id}"), None).await;
        if record["status"] == json!(status) {
          return record;
        }
        tokio::time::sleep(Duration::from_millis(100)).await;
      }
      panic!("run {id} did not reach status {status}");
    }

    pub(super) async fn read_events(response: Response) -> Vec<SseEvent> {
      let mut stream = response.into_body().into_data_stream();
      let mut buffer = String::new();
      let mut events = vec![];
      while let Some(event) = next_event(&mut stream, &mut buffer).await {
        events.push(event);
      }
      events
    }

    pub(super) async fn next_event(stream: &mut BodyDataStream, buffer: &mut String) -> Option<SseEvent> {
      loop {
        if let Some(end) = buffer.find("\n\n") {
          let block: String = buffer.drain(..end + 2).collect();
          if let Some(event) = parse_event(&block) {
            return Some(event);
          }
          continue;
        }
        let chunk = stream.next().await?.unwrap();
        buffer.push_str(str::from_utf8(&chunk).unwrap());
      }
    }

    fn parse_event(block: &str) -> Option<SseEvent> {
      let mut name = None;
      let mut data = String::new();
      for line in block.lines() {
        if let Some(value) = line.strip_prefix("event: ") {
          name = Some(value.to_owned());
        } else if let Some(value) = line.strip_prefix("data: ") {
          data.push_str(value);
        }
      }
      Some(SseEvent {
        name: name?,
        data: serde_json::from_str(&data).unwrap(),
      })
    }
  }
}
