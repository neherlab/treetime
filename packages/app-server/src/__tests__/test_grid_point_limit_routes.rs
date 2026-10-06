#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{
    app_with_grid_limit, events_of, request, timetree_config, wait_for_status,
  };
  use helpers::with_limit;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_grid::MaxGridPoints;

  const ABOVE_SERVER_LIMIT: &str = "`max_grid_points` is 3000, more than the limit of this server, 2000 \
    (treetime-server --max-grid-points); set 2000 or less, or leave the setting unset";

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_grid_point_limit_routes_preview_fills_the_server_limit() {
    let test = app_with_grid_limit(MaxGridPoints::new(2_000).unwrap());
    let (_, response) = request(
      &test,
      "POST",
      "/api/run-config",
      Some(json!({ "command": "timetree", "config": timetree_config() })),
    )
    .await;
    let flag = response["code"]["command_line_text"]
      .as_str()
      .is_some_and(|text| text.contains("--max-grid-points 2000"));
    assert_eq!(
      (json!("valid"), json!(2000), true),
      (
        response["status"].clone(),
        response["config"]["max_grid_points"].clone(),
        flag
      ),
      "{response:#}"
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_grid_point_limit_routes_preview_keeps_a_lower_value() {
    let test = app_with_grid_limit(MaxGridPoints::new(2_000).unwrap());
    let (_, response) = request(
      &test,
      "POST",
      "/api/run-config",
      Some(json!({ "command": "timetree", "config": with_limit(1_500) })),
    )
    .await;
    assert_eq!(
      (json!("valid"), json!(1500)),
      (
        response["status"].clone(),
        response["config"]["max_grid_points"].clone()
      )
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_grid_point_limit_routes_preview_rejects_a_value_above_the_server_limit() {
    let test = app_with_grid_limit(MaxGridPoints::new(2_000).unwrap());
    let (_, response) = request(
      &test,
      "POST",
      "/api/run-config",
      Some(json!({ "command": "timetree", "config": with_limit(3_000) })),
    )
    .await;
    assert_eq!(
      (json!("invalid"), json!(ABOVE_SERVER_LIMIT)),
      (response["status"].clone(), response["message"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_grid_point_limit_routes_run_above_the_server_limit_ends_with_an_error() {
    let test = app_with_grid_limit(MaxGridPoints::new(2_000).unwrap());
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "timetree", "config": with_limit(3_000) })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let events = events_of(&test, &id, "").await;
    let terminal = &events.last().unwrap().data["data"];
    assert_eq!(
      (json!("error"), json!(ABOVE_SERVER_LIMIT)),
      (terminal["status"].clone(), terminal["message"].clone())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_grid_point_limit_routes_run_records_the_filled_limit_as_a_changed_setting() {
    let test = app_with_grid_limit(MaxGridPoints::default());
    let (_, record) = request(
      &test,
      "POST",
      "/api/runs",
      Some(json!({ "command": "timetree", "config": timetree_config() })),
    )
    .await;
    let id = record["id"].as_str().unwrap().to_owned();
    let record = wait_for_status(&test, &id, "ok").await;
    let changed = record["changed_settings"]
      .as_array()
      .is_some_and(|settings| settings.contains(&json!("max_grid_points")));
    assert_eq!(
      (json!(1_000_000), true),
      (record["config"]["max_grid_points"].clone(), changed),
      "{record:#}"
    );
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::timetree_config;
    use serde_json::{Value, json};

    pub(super) fn with_limit(max_grid_points: usize) -> Value {
      let mut config = timetree_config();
      config["max_grid_points"] = json!(max_grid_points);
      config
    }
  }
}
