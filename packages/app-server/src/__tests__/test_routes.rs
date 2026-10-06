#[cfg(test)]
pub(crate) mod tests {
  use helpers::{
    TestApp, app, app_with_upload_limit, create_deferred, events_of, next_event, read_events, request, timetree_config,
    wait_for_status,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_run_streams_events_and_records_outputs() {
    let test = app();
    let (status, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "timetree", "config": timetree_config() })),
    )
    .await;
    assert_eq!((200, json!("running")), (status, record["status"].clone()));
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
          "nwk",
          "nexus",
          "auspice",
          "augur-node-data",
          "gtr",
          "reconstructed-nuc-fasta",
          "clock-model",
          "coalescent-tsv",
          "tracelog",
          "clock-csv"
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

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown_key(  json!({ "tree": "zika/20/tree.nwk", "bogus": 1 }),      "invalid configuration: unknown field `bogus`")]
  #[case::null_value(   json!({ "tree": null }),                                 "invalid configuration: null is not of type \"string\"")]
  #[case::wrong_type(   json!({ "tree": 7 }),                                    "invalid configuration: 7 is not of type \"string\"")]
  #[trace]
  #[tokio::test]
  async fn test_routes_rejected_configuration_answers_an_invalid_request(
    #[case] config: Value,
    #[case] message: &str,
  ) {
    let test = app();
    let (status, error) = request(&test, "POST", "/api/runs", Some(json!({ "command": "clock", "config": config }))).await;
    let (_, list) = request(&test, "GET", "/api/runs", None).await;
    assert_eq!(
      (400, json!("invalid_request"), json!(message), json!([])),
      (status, error["code"].clone(), error["message"].clone(), list["runs"].clone())
    );
  }

  #[tokio::test]
  async fn test_routes_start_rejects_a_configuration_that_does_not_parse() {
    let test = app();
    let id = create_deferred(&test).await;
    let (status, error) = request(
      &test,
      "POST",
      &format!("/api/runs/{id}/start"),
      Some(json!({ "config": { "tree": "zika/20/tree.nwk", "bogus": 1 } })),
    )
    .await;
    let (_, record) = request(&test, "GET", &format!("/api/runs/{id}"), None).await;
    assert_eq!(
      (
        400,
        json!("invalid configuration: unknown field `bogus`"),
        json!("created")
      ),
      (status, error["message"].clone(), record["status"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_input_outside_the_examples_dir_ends_with_an_error_run() {
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
  async fn test_routes_create_run_takes_no_title() {
    let test = app();
    let (status, _) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "clock", "config": {}, "title": "zika", "defer_start": true })),
    )
    .await;
    assert_eq!(
      400, status,
      "a run gets its title from the server and is renamed afterwards"
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_list_rename_and_pin() {
    let test = app();
    let ids = [create_deferred(&test).await, create_deferred(&test).await];
    let (_, list) = request(&test, "GET", "/api/runs", None).await;
    let listed: Vec<&str> = list["runs"]
      .as_array()
      .unwrap()
      .iter()
      .map(|run| run["id"].as_str().unwrap())
      .collect();
    assert_eq!(
      (vec![ids[1].as_str(), ids[0].as_str()], json!(0)),
      (listed, list["active_runs"].clone()),
      "the list starts with the newest run"
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
    let data = TestApp::examples_dir();
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
        "command": "timetree",
        "config": {
          "tree": "zika/86/tree.nwk",
          "metadata": "zika/86/metadata.tsv",
          "alignment": ["zika/86/aln.fasta.xz"],
        },
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
  async fn test_routes_check_inputs_rejects_paths_outside_the_examples_dir() {
    let test = app();
    let (status, body) = request(
      &test,
      "POST",
      "/api/check-inputs",
      Some(json!({ "command": "prune", "config": { "tree": "../Cargo.toml" } })),
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
        json!("invalid configuration: \"fast\" is not of type \"number\"")
      ),
      (response["status"].clone(), response["message"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_run_config_adds_the_run_outputs() {
    let test = app();
    let (_, response) = request(
      &test,
      "POST",
      "/api/run-config",
      Some(
        json!({ "command": "clock", "config": { "tree": "t.nwk", "metadata": "m.tsv", "output_selection": ["nwk"] } }),
      ),
    )
    .await;
    assert_eq!(
      (json!("valid"), json!(["nwk", "auspice"]), json!("out")),
      (
        response["status"].clone(),
        response["config"]["output_selection"].clone(),
        response["config"]["output_all"].clone()
      )
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
  async fn test_routes_results_comparison_and_clades_of_a_finished_run() {
    let test = app();
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "timetree", "config": timetree_config() })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let (early, _) = request(&test, "GET", &format!("/api/runs/{id}/results"), None).await;
    wait_for_status(&test, &id, "ok").await;

    let (status, results) = request(&test, "GET", &format!("/api/runs/{id}/results"), None).await;
    let (_, comparison) = request(&test, "GET", &format!("/api/runs/{id}/compare/{id}"), None).await;
    let root = results["tree"]["nodes"][0]["name"].clone();
    let (_, clade) = request(
      &test,
      "POST",
      "/api/clade-in-runs",
      Some(json!({ "run": id, "node": root })),
    )
    .await;

    assert_eq!(
      (
        true,
        200,
        json!("timetree"),
        json!(20),
        json!(0.0),
        json!({ "matches": [], "searched_runs": 0, "unreadable_runs": [] })
      ),
      (
        early == 200 || early == 409,
        status,
        results["results"]["command"].clone(),
        results["results"]["data"]["estimates"]["samples"].clone(),
        comparison["estimates"]["root_shift_days"].clone(),
        clade
      )
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_invalid_run_id_is_a_bad_request() {
    let test = app();
    let (status, _) = request(&test, "GET", "/api/runs/a.b", None).await;
    assert_eq!(400, status);
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_malformed_json_body_is_a_bad_request() {
    let test = app();
    let response = test
      .send(
        axum::http::Request::post("/api/runs")
          .header("content-type", "application/json")
          .body(axum::body::Body::from("{\"command\":"))
          .unwrap(),
      )
      .await;
    assert_eq!(
      (
        400,
        json!({
          "code": "invalid_request",
          "message": "Failed to parse the request body as JSON: EndOfFile: unexpected end of file at line 1 column 12",
          "causes": [],
        })
      ),
      helpers::body_json(response).await
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_malformed_query_parameter_is_a_bad_request() {
    let test = app();
    let (status, body) = request(&test, "GET", "/api/runs/abc/events?from=x", None).await;
    assert_eq!(
      (
        400,
        json!({
          "code": "invalid_request",
          "message": "Failed to deserialize query string: Unexpected: invalid value \"x\", expected usize at line 1 column 6 (path: from)",
          "causes": [],
        })
      ),
      (status, body)
    );
  }

  pub(crate) mod helpers {
    use crate::create_router;
    use crate::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig};
    use crate::web::WebOptions;
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
    use tokio_util::sync::CancellationToken;
    use tower::ServiceExt;
    use treetime_grid::MaxGridPoints;

    pub(crate) struct SseEvent {
      pub id: Option<String>,
      pub name: String,
      pub data: Value,
    }

    pub(crate) struct TestApp {
      pub router: Router,
      pub runs_dir: TempDir,
    }

    impl TestApp {
      pub(crate) fn examples_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data")
      }

      pub(crate) async fn send(&self, request: Request<Body>) -> Response {
        self.router.clone().oneshot(request).await.unwrap()
      }

      pub(crate) async fn upload(&self, id: &str, name: &str, bytes: Vec<u8>) -> (u16, Value) {
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

    pub(crate) fn app() -> TestApp {
      app_with_upload_limit(DEFAULT_MAX_UPLOAD_SIZE)
    }

    pub(crate) fn app_with_upload_limit(max_upload_size: usize) -> TestApp {
      app_with(max_upload_size, &CancellationToken::new(), &WebOptions::default())
    }

    pub(crate) fn app_with(max_upload_size: usize, shutdown: &CancellationToken, options: &WebOptions) -> TestApp {
      app_from(max_upload_size, None, shutdown, options)
    }

    pub(crate) fn app_with_grid_limit(max_grid_points: MaxGridPoints) -> TestApp {
      app_from(
        DEFAULT_MAX_UPLOAD_SIZE,
        Some(max_grid_points),
        &CancellationToken::new(),
        &WebOptions::default(),
      )
    }

    fn app_from(
      max_upload_size: usize,
      max_grid_points: Option<MaxGridPoints>,
      shutdown: &CancellationToken,
      options: &WebOptions,
    ) -> TestApp {
      let runs_dir = tempdir().unwrap();
      let router = create_router(
        ServerConfig {
          examples_dir: TestApp::examples_dir(),
          runs_dir: runs_dir.path().to_path_buf(),
          max_upload_size,
          max_grid_points,
          shutdown: shutdown.clone(),
          settings: None,
        },
        options,
      )
      .unwrap();
      TestApp { router, runs_dir }
    }

    pub(crate) fn timetree_config() -> Value {
      json!({
        "tree": "zika/20/tree.nwk",
        "metadata": "zika/20/metadata.tsv",
        "alignment": ["zika/20/aln.fasta.xz"],
        "max_iter": 2,
        "seed": 7,
      })
    }

    pub(crate) async fn request(test: &TestApp, method: &str, uri: &str, body: Option<Value>) -> (u16, Value) {
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

    pub(crate) async fn body_json(response: Response) -> (u16, Value) {
      let status = response.status().as_u16();
      let bytes = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
      let value = if bytes.is_empty() {
        Value::Null
      } else {
        serde_json::from_slice(&bytes).unwrap_or_else(|_| json!(String::from_utf8_lossy(&bytes)))
      };
      (status, value)
    }

    pub(crate) async fn events_of(test: &TestApp, id: &str, query: &str) -> Vec<SseEvent> {
      let response = test
        .send(
          Request::get(format!("/api/runs/{id}/events{query}"))
            .body(Body::empty())
            .unwrap(),
        )
        .await;
      read_events(response).await
    }

    pub(crate) async fn open_app_events(test: &TestApp, query: &str, last_event_id: Option<&str>) -> BodyDataStream {
      let request = Request::get(format!("/api/events{query}"));
      let request = match last_event_id {
        Some(id) => request.header("last-event-id", id),
        None => request,
      };
      let response = test.send(request.body(Body::empty()).unwrap()).await;
      assert_eq!(
        (200, Some("text/event-stream")),
        (
          response.status().as_u16(),
          response
            .headers()
            .get("content-type")
            .and_then(|value| value.to_str().ok())
        )
      );
      response.into_body().into_data_stream()
    }

    pub(crate) async fn take_events(stream: &mut BodyDataStream, count: usize) -> Vec<SseEvent> {
      let mut buffer = String::new();
      let mut events = vec![];
      while events.len() < count {
        let event = tokio::time::timeout(Duration::from_secs(10), next_event(stream, &mut buffer))
          .await
          .expect("an app event within 10 s")
          .expect("an open app event stream");
        events.push(event);
      }
      events
    }

    pub(crate) async fn create_deferred(test: &TestApp) -> String {
      let (_, record) = request(
        test,
        "POST",
        "/api/runs",
        Some(json!({ "command": "clock", "config": {}, "defer_start": true })),
      )
      .await;
      record["id"].as_str().unwrap().to_owned()
    }

    pub(crate) fn is_documented_path(doc: &Value, path: &str) -> bool {
      is_documented(doc, path, |template, segments| template == segments)
    }

    pub(crate) fn is_documented_prefix(doc: &Value, path: &str) -> bool {
      is_documented(doc, path, |template, segments| template >= segments)
    }

    fn is_documented(doc: &Value, path: &str, lengths_match: impl Fn(usize, usize) -> bool) -> bool {
      let segments = path.split('/').collect::<Vec<_>>();
      doc["paths"].as_object().unwrap().keys().any(|template| {
        let template = template.split('/').collect::<Vec<_>>();
        lengths_match(template.len(), segments.len())
          && segments
            .iter()
            .zip(&template)
            .all(|(segment, pattern)| segment == pattern || pattern.starts_with('{'))
      })
    }

    pub(crate) async fn wait_for_status(test: &TestApp, id: &str, status: &str) -> Value {
      for _ in 0..600 {
        let (_, record) = request(test, "GET", &format!("/api/runs/{id}"), None).await;
        if record["status"] == json!(status) {
          return record;
        }
        tokio::time::sleep(Duration::from_millis(100)).await;
      }
      panic!("run {id} did not reach status {status}");
    }

    pub(crate) async fn read_events(response: Response) -> Vec<SseEvent> {
      let mut stream = response.into_body().into_data_stream();
      let mut buffer = String::new();
      let mut events = vec![];
      while let Some(event) = next_event(&mut stream, &mut buffer).await {
        events.push(event);
      }
      events
    }

    pub(crate) async fn next_event(stream: &mut BodyDataStream, buffer: &mut String) -> Option<SseEvent> {
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
      let mut id = None;
      let mut name = None;
      let mut data = String::new();
      for line in block.lines() {
        if let Some(value) = line.strip_prefix("id: ") {
          id = Some(value.to_owned());
        } else if let Some(value) = line.strip_prefix("event: ") {
          name = Some(value.to_owned());
        } else if let Some(value) = line.strip_prefix("data: ") {
          data.push_str(value);
        }
      }
      Some(SseEvent {
        id,
        name: name?,
        data: serde_json::from_str(&data).unwrap(),
      })
    }
  }
}
