#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{
    SseEvent, app, create_deferred, is_documented_path, is_documented_prefix, open_app_events, request, take_events,
  };
  use crate::routes::api_doc;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};

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
    let seqs = |events: &[SseEvent]| events.iter().map(|event| event.data["seq"].clone()).collect::<Vec<_>>();
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
            { "path": "/api/datasets", "scope": "exact" },
            { "path": "/api/examples/download", "scope": "exact" },
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
}
