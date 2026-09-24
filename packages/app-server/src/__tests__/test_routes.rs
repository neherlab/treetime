#[cfg(test)]
mod tests {
  use crate::routes::api_doc;
  use helpers::{SseEvent, app, check_config, next_event, post_json, read_events, timetree_config};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use tempfile::tempdir;

  #[tokio::test]
  async fn test_routes_command_streams_started_progress_and_one_terminal() {
    let out = tempdir().unwrap();
    let app = app(out.path());
    let events = read_events(post_json(&app, "/api/timetree", &timetree_config()).await).await;

    let names: Vec<&str> = events.iter().map(|event| event.name.as_str()).collect();
    assert_eq!(Some(&"started"), names.first());
    assert_eq!(Some(&"terminal"), names.last());
    assert_eq!(1, names.iter().filter(|name| **name == "terminal").count());
    assert!(names.contains(&"progress"));

    let job_id = events[0].data["job_id"].as_str().unwrap().to_owned();
    let terminal = &events[events.len() - 1].data;
    assert_eq!(
      (json!("ok"), json!(job_id), json!("timetree")),
      (
        terminal["status"].clone(),
        terminal["job_id"].clone(),
        terminal["result"]["command"].clone()
      )
    );
    let outputs = terminal["result"]["output_files"].as_array().unwrap();
    assert!(
      outputs.iter().all(|path| path
        .as_str()
        .unwrap()
        .starts_with(&out.path().join(&job_id).to_string_lossy().into_owned())),
      "outputs go to the job's own directory: {outputs:?}"
    );
  }

  #[tokio::test]
  async fn test_routes_unknown_setting_ends_with_error_terminal() {
    let out = tempdir().unwrap();
    let app = app(out.path());
    let events =
      read_events(post_json(&app, "/api/clock", &json!({ "tree": "zika/20/tree.nwk", "bogus": 1 })).await).await;
    let terminal = &events[events.len() - 1];
    assert_eq!(
      (
        2,
        "terminal",
        json!("error"),
        json!("invalid configuration: unknown field `bogus`")
      ),
      (
        events.len(),
        terminal.name.as_str(),
        terminal.data["status"].clone(),
        terminal.data["message"].clone()
      )
    );
  }

  #[tokio::test]
  async fn test_routes_input_outside_data_dir_ends_with_error_terminal() {
    let out = tempdir().unwrap();
    let app = app(out.path());
    let config = json!({ "tree": "../Cargo.toml", "metadata": "zika/20/metadata.tsv" });
    let events = read_events(post_json(&app, "/api/clock", &config).await).await;
    let terminal = &events[events.len() - 1].data;
    assert_eq!(
      (
        json!("error"),
        json!("input `../Cargo.toml` of setting `tree` is outside the directories the server reads inputs from")
      ),
      (terminal["status"].clone(), terminal["message"].clone())
    );
  }

  #[tokio::test]
  async fn test_routes_cancel_one_of_two_concurrent_jobs() {
    let out = tempdir().unwrap();
    let app = app(out.path());

    let mut first = post_json(&app, "/api/timetree", &timetree_config())
      .await
      .into_body()
      .into_data_stream();
    let mut second = post_json(&app, "/api/timetree", &timetree_config())
      .await
      .into_body()
      .into_data_stream();
    let mut first_buffer = String::new();
    let started = next_event(&mut first, &mut first_buffer).await.unwrap();
    let first_id = started.data["job_id"].as_str().unwrap().to_owned();

    let response = post_json(&app, &format!("/api/jobs/{first_id}/cancel"), &json!({})).await;
    assert_eq!(200, response.status().as_u16());

    let mut first_events = vec![started];
    while let Some(event) = next_event(&mut first, &mut first_buffer).await {
      first_events.push(event);
    }
    let mut second_buffer = String::new();
    let mut second_events: Vec<SseEvent> = vec![];
    while let Some(event) = next_event(&mut second, &mut second_buffer).await {
      second_events.push(event);
    }

    let terminals = |events: &[SseEvent]| -> Vec<String> {
      events
        .iter()
        .filter(|event| event.name == "terminal")
        .map(|event| event.data["status"].as_str().unwrap().to_owned())
        .collect()
    };
    assert_eq!(
      (vec!["cancelled".to_owned()], vec!["ok".to_owned()]),
      (terminals(&first_events), terminals(&second_events))
    );

    let response = post_json(&app, &format!("/api/jobs/{first_id}/cancel"), &json!({})).await;
    assert_eq!(404, response.status().as_u16(), "a finished job cannot be cancelled");
  }

  #[tokio::test]
  async fn test_routes_cancel_rejects_invalid_job_id() {
    let out = tempdir().unwrap();
    let app = app(out.path());
    let response = post_json(&app, "/api/jobs/a.b/cancel", &json!({})).await;
    assert_eq!(400, response.status().as_u16());
  }

  #[tokio::test]
  async fn test_routes_check_config_reports_cli_errors() {
    let out = tempdir().unwrap();
    let app = app(out.path());
    let response = check_config(
      &app,
      json!({ "command": "prune", "text": "tree: t.nwk\nprune_short: fast\n" }),
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

  #[tokio::test]
  async fn test_routes_check_config_returns_normalized_config() {
    let out = tempdir().unwrap();
    let app = app(out.path());
    let response = check_config(&app, json!({ "command": "timetree", "text": "tree: t.nwk\n" })).await;
    assert_eq!(
      (json!("valid"), json!("t.nwk"), json!(3.0)),
      (
        response["status"].clone(),
        response["config"]["tree"].clone(),
        response["config"]["clock_filter"].clone()
      )
    );
  }

  #[test]
  fn test_routes_openapi_describes_config_bodies_and_bridge_types() {
    let doc = api_doc().unwrap();
    let operation = &doc["paths"]["/api/timetree"]["post"];
    assert_eq!(
      (
        json!("command"),
        json!("#/components/schemas/TimetreeConfig"),
        json!("#/components/schemas/JobEvent"),
        json!("request"),
        json!("#/components/schemas/CheckConfigRequest"),
        json!("query"),
      ),
      (
        operation["x-bridge-type"].clone(),
        operation["requestBody"]["content"]["application/json"]["schema"]["$ref"].clone(),
        operation["responses"]["200"]["content"]["text/event-stream"]["schema"]["$ref"].clone(),
        doc["paths"]["/api/check-config"]["post"]["x-bridge-type"].clone(),
        doc["paths"]["/api/check-config"]["post"]["requestBody"]["content"]["application/json"]["schema"]["$ref"]
          .clone(),
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
    use crate::state::ServerConfig;
    use axum::Router;
    use axum::body::{Body, BodyDataStream};
    use axum::http::Request;
    use axum::response::Response;
    use serde_json::{Value, json};
    use std::path::{Path, PathBuf};
    use std::str;
    use tokio_stream::StreamExt;
    use tower::ServiceExt;

    pub(super) struct SseEvent {
      pub name: String,
      pub data: Value,
    }

    pub(super) fn app(out_dir: &Path) -> Router {
      let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data");
      create_router(
        ServerConfig {
          data_dir,
          out_dir: out_dir.to_path_buf(),
        },
        None,
      )
      .unwrap()
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

    pub(super) async fn post_json(app: &Router, uri: &str, body: &Value) -> Response {
      let request = Request::post(uri)
        .header("content-type", "application/json")
        .body(Body::from(body.to_string()))
        .unwrap();
      app.clone().oneshot(request).await.unwrap()
    }

    pub(super) async fn check_config(app: &Router, body: Value) -> Value {
      let response = post_json(app, "/api/check-config", &body).await;
      let bytes = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
      serde_json::from_slice(&bytes).unwrap()
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
