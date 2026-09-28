#[cfg(test)]
mod tests {
  use crate::routes::api_doc;
  use helpers::{
    TestApp, app, app_with_upload_limit, create_deferred, events_of, is_documented_path, is_documented_prefix,
    next_event, open_app_events, read_events, request, take_events, timetree_config, wait_for_status,
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
          "Nwk",
          "Nexus",
          "Auspice",
          "AugurNodeData",
          "Gtr",
          "ReconstructedNucFasta",
          "ClockModel",
          "CoalescentTsv",
          "Tracelog",
          "ClockCsv"
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
  async fn test_routes_check_inputs_rejects_paths_outside_the_data_dir() {
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
        json!("invalid configuration: \"fast\" is not of types \"null\", \"number\"")
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
        json!({ "command": "clock", "config": { "tree": "t.nwk", "metadata": "m.tsv", "output_selection": ["Nwk"] } }),
      ),
    )
    .await;
    assert_eq!(
      (json!("valid"), json!(["Nwk", "Auspice"]), json!("out")),
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
          "message": "Failed to parse the request body as JSON: command: EOF while parsing a value at line 1 column 11",
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
          "message": "Failed to deserialize query string: from: invalid digit found in string",
          "causes": [],
        })
      ),
      (status, body)
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_app_events_report_run_changes_with_their_stale_paths() {
    let test = app();
    let mut stream = open_app_events(&test, "", None).await;
    let id = create_deferred(&test).await;
    request(
      &test,
      "PATCH",
      &format!("/api/runs/{id}"),
      Some(json!({ "title": "renamed" })),
    )
    .await;
    request(
      &test,
      "PATCH",
      &format!("/api/runs/{id}"),
      Some(json!({ "pinned": true })),
    )
    .await;

    let events = take_events(&mut stream, 3).await;
    let stale = json!([
      { "path": "/api/runs", "scope": "exact" },
      { "path": format!("/api/runs/{id}"), "scope": "subtree" },
      { "path": "/api/clade-in-runs", "scope": "exact" },
    ]);
    assert_eq!(
      vec![
        ("run-created", json!("run-created"), json!(id), stale.clone()),
        ("run-updated", json!("run-updated"), json!(id), stale.clone()),
        ("run-updated", json!("run-updated"), json!(id), stale),
      ],
      events
        .iter()
        .map(|event| (
          event.name.as_str(),
          event.data["kind"].clone(),
          event.data["run"]["id"].clone(),
          event.data["stale"].clone(),
        ))
        .collect::<Vec<_>>()
    );
    assert_eq!(
      (json!("renamed"), json!(true)),
      (
        events[1].data["run"]["title"].clone(),
        events[2].data["run"]["pinned"].clone()
      )
    );
    let seqs = events
      .iter()
      .map(|event| event.data["seq"].as_u64().unwrap())
      .collect::<Vec<_>>();
    assert_eq!(
      (
        vec![seqs[0], seqs[0] + 1, seqs[0] + 2],
        seqs.iter().map(u64::to_string).collect::<Vec<_>>()
      ),
      (
        seqs.clone(),
        events.iter().map(|event| event.id.clone().unwrap()).collect::<Vec<_>>()
      ),
      "consecutive numbers, sent as the SSE event id"
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_app_events_resume_from_an_event_by_query_or_last_event_id() {
    let test = app();
    let mut live = open_app_events(&test, "", None).await;
    for _ in 0..3 {
      create_deferred(&test).await;
    }
    let seqs = |events: &[helpers::SseEvent]| events.iter().map(|event| event.data["seq"].clone()).collect::<Vec<_>>();
    let all = take_events(&mut live, 3).await;
    let first = all[0].data["seq"].as_u64().unwrap();

    let mut from = open_app_events(&test, &format!("?from={}", first + 1), None).await;
    let mut last_event_id = open_app_events(&test, "", Some(&first.to_string())).await;
    assert_eq!(
      (seqs(&all[1..]), seqs(&all[1..])),
      (
        seqs(&take_events(&mut from, 2).await),
        seqs(&take_events(&mut last_event_id, 2).await)
      )
    );

    create_deferred(&test).await;
    let next = json!(first + 3);
    assert_eq!(
      (vec![next.clone()], vec![next]),
      (
        seqs(&take_events(&mut from, 1).await),
        seqs(&take_events(&mut last_event_id, 1).await)
      ),
      "a resumed stream follows new changes"
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_app_events_resync_when_the_event_is_not_in_the_log() {
    let test = app();
    let mut live = open_app_events(&test, "", None).await;
    create_deferred(&test).await;
    let mut head = take_events(&mut live, 1).await[0].data["seq"].as_u64().unwrap();

    for query in ["?from=1".to_owned(), format!("?from={}", head + 1000)] {
      let mut stream = open_app_events(&test, &query, None).await;
      let resync = &take_events(&mut stream, 1).await[0];
      assert_eq!(
        (
          "resync",
          Some(head.to_string()),
          json!(head),
          json!([
            { "path": "/api/runs", "scope": "subtree" },
            { "path": "/api/clade-in-runs", "scope": "exact" },
          ])
        ),
        (
          resync.name.as_str(),
          resync.id.clone(),
          resync.data["seq"].clone(),
          resync.data["stale"].clone()
        ),
        "{query}"
      );
      create_deferred(&test).await;
      let next = &take_events(&mut stream, 1).await[0];
      let expected = take_events(&mut live, 1).await[0].data["seq"].clone();
      assert_eq!(
        (json!("run-created"), expected.clone()),
        (next.data["kind"].clone(), next.data["seq"].clone())
      );
      head = expected.as_u64().unwrap();
    }
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_routes_app_event_stale_paths_are_paths_or_path_prefixes_of_the_api() {
    let test = app();
    let mut stream = open_app_events(&test, "", None).await;
    let id = create_deferred(&test).await;
    request(
      &test,
      "PATCH",
      &format!("/api/runs/{id}"),
      Some(json!({ "title": "renamed" })),
    )
    .await;
    let mut stale = take_events(&mut stream, 2)
      .await
      .iter()
      .flat_map(|event| event.data["stale"].as_array().unwrap().clone())
      .collect::<Vec<_>>();
    let mut resync = open_app_events(&test, "?from=1", None).await;
    stale.extend(
      take_events(&mut resync, 1).await[0].data["stale"]
        .as_array()
        .unwrap()
        .clone(),
    );
    let doc = api_doc().unwrap();
    let undocumented = stale
      .iter()
      .filter(|entry| {
        let path = entry["path"].as_str().unwrap();
        match entry["scope"].as_str().unwrap() {
          "exact" => !is_documented_path(&doc, path),
          _ => !is_documented_prefix(&doc, path),
        }
      })
      .collect::<Vec<_>>();
    assert_eq!((false, Vec::<&Value>::new()), (stale.is_empty(), undocumented));
  }

  #[test]
  fn test_routes_openapi_documents_the_app_event_stream() {
    let doc = api_doc().unwrap();
    let operation = &doc["paths"]["/api/events"]["get"];
    assert_eq!(
      (
        json!("events"),
        json!("#/components/schemas/AppEvent"),
        json!(["from"]),
        json!("kind"),
      ),
      (
        operation["operationId"].clone(),
        operation["responses"]["200"]["content"]["text/event-stream"]["schema"]["$ref"].clone(),
        json!(
          operation["parameters"]
            .as_array()
            .unwrap()
            .iter()
            .map(|parameter| parameter["name"].clone())
            .collect::<Vec<_>>()
        ),
        doc["components"]["schemas"]["AppEvent"]["discriminator"]["propertyName"].clone(),
      )
    );
  }

  #[test]
  fn test_routes_openapi_operations_have_ids_descriptions_and_error_responses() {
    let doc = api_doc().unwrap();
    let incomplete = doc["paths"]
      .as_object()
      .unwrap()
      .iter()
      .flat_map(|(path, item)| {
        item
          .as_object()
          .unwrap()
          .iter()
          .map(move |(method, operation)| (format!("{method} {path}"), operation))
      })
      .filter(|(_, operation)| {
        !operation["operationId"].is_string()
          || operation["description"].as_str().is_none_or(str::is_empty)
          || operation["responses"]["default"]["content"]["application/json"]["schema"]["$ref"]
            != json!("#/components/schemas/ErrorResponse")
      })
      .map(|(operation, _)| operation)
      .collect::<Vec<_>>();
    assert_eq!(Vec::<String>::new(), incomplete);
  }

  #[test]
  fn test_routes_openapi_describes_operations_and_their_bodies() {
    let doc = api_doc().unwrap();
    assert_eq!(
      (
        json!("runsCreate"),
        json!("#/components/schemas/CreateRunRequest"),
        json!("#/components/schemas/ErrorResponse"),
        json!("#/components/schemas/RunEvent"),
      ),
      (
        doc["paths"]["/api/runs"]["post"]["operationId"].clone(),
        doc["paths"]["/api/runs"]["post"]["requestBody"]["content"]["application/json"]["schema"]["$ref"].clone(),
        doc["paths"]["/api/runs/{id}"]["get"]["responses"]["default"]["content"]["application/json"]["schema"]["$ref"]
          .clone(),
        doc["paths"]["/api/runs/{id}/events"]["get"]["responses"]["200"]["content"]["text/event-stream"]["schema"]
          ["$ref"]
          .clone(),
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

  #[test]
  fn test_routes_openapi_config_components_carry_the_cli_annotations() {
    let doc = api_doc().unwrap();
    let schemas = &doc["components"]["schemas"];
    assert_eq!(
      (json!("--relax"), json!("--branch-split-grid-n-points"), json!("input"),),
      (
        schemas["TimetreeConfig"]["properties"]["relax"]["x-cli-flag"].clone(),
        schemas["BranchSplitArgs"]["properties"]["n_points"]["x-cli-flag"].clone(),
        schemas["ClockConfig"]["properties"]["tree"]["x-path"].clone(),
      )
    );
  }

  #[test]
  fn test_routes_openapi_carries_the_setting_catalog_of_every_command() {
    let doc = api_doc().unwrap();
    let commands = doc["x-setting-catalog"]["commands"]
      .as_array()
      .unwrap()
      .iter()
      .map(|settings| settings["command"].clone())
      .collect::<Vec<_>>();
    assert_eq!(
      (
        vec![
          json!("timetree"),
          json!("optimize"),
          json!("prune"),
          json!("ancestral"),
          json!("clock"),
          json!("mugration")
        ],
        true
      ),
      (commands, doc["components"]["schemas"]["SettingCatalog"].is_object())
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
      pub id: Option<String>,
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

    pub(super) async fn open_app_events(test: &TestApp, query: &str, last_event_id: Option<&str>) -> BodyDataStream {
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

    pub(super) async fn take_events(stream: &mut BodyDataStream, count: usize) -> Vec<SseEvent> {
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

    pub(super) async fn create_deferred(test: &TestApp) -> String {
      let (_, record) = request(
        test,
        "POST",
        "/api/runs",
        Some(json!({ "command": "clock", "config": {}, "defer_start": true })),
      )
      .await;
      record["id"].as_str().unwrap().to_owned()
    }

    pub(super) fn is_documented_path(doc: &Value, path: &str) -> bool {
      is_documented(doc, path, |template, segments| template == segments)
    }

    pub(super) fn is_documented_prefix(doc: &Value, path: &str) -> bool {
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
