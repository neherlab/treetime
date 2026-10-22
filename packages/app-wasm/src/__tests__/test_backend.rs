#[cfg(test)]
mod tests {
  use crate::backend::{header_list, web_router};
  use crate::exchange::PortExchange;
  use app_backend::bridge::executor::ThreadExecutor;
  use app_server::port::port_request;
  use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig};
  use helpers::{config, fetch_json};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::path::PathBuf;
  use std::sync::Arc;
  use tempfile::tempdir;
  use tokio_util::sync::CancellationToken;
  use treetime_schema::version_info;
  use treetime_utils::io::json::to_json_value;
  use treetime_utils::{assert_error, o, vec_of_owned};

  #[tokio::test]
  async fn test_backend_router_answers_the_public_api() {
    let runs = tempdir().unwrap();
    let router = web_router(config(runs.path().to_owned()), Arc::new(ThreadExecutor)).unwrap();
    let expected = (200, to_json_value(&version_info()).unwrap());
    assert_eq!(
      expected,
      fetch_json(PortExchange::start(
        router,
        port_request("GET", "/api/version", vec![], None)
      ))
      .await
    );
  }

  #[tokio::test]
  async fn test_backend_router_has_no_routes_of_a_local_app() {
    let runs = tempdir().unwrap();
    let router = web_router(config(runs.path().to_owned()), Arc::new(ThreadExecutor)).unwrap();
    let (status, body) = fetch_json(PortExchange::start(
      router,
      port_request("GET", "/api/app-settings", vec![], None),
    ))
    .await;
    assert_eq!((404, json!("not_found")), (status, body["code"].clone()));
  }

  #[test]
  fn test_backend_header_list_pairs_names_with_values() {
    let expected = vec![(o!("accept"), o!("text/event-stream")), (o!("x-a"), o!("1"))];
    let actual = header_list(vec_of_owned!["accept", "x-a"], vec_of_owned!["text/event-stream", "1"]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_backend_header_list_rejects_unequal_lengths() {
    assert_error!(
      header_list(vec_of_owned!["accept", "x-a"], vec_of_owned!["text/event-stream"]),
      "a request has 2 header names and 1 header values"
    );
  }

  mod helpers {
    use super::{CancellationToken, DEFAULT_MAX_UPLOAD_SIZE, PathBuf, PortExchange, ServerConfig, Value};

    pub(super) fn config(runs_dir: PathBuf) -> ServerConfig {
      ServerConfig {
        runs_dir,
        max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
        max_grid_points: None,
        shutdown: CancellationToken::new(),
        settings: None,
      }
    }

    pub(super) async fn fetch_json(exchange: PortExchange) -> (u16, Value) {
      let head = exchange.head().await.unwrap();
      let mut body = vec![];
      while let Some(chunk) = exchange.chunk().await.unwrap() {
        body.extend_from_slice(&chunk);
      }
      (head.status, serde_json::from_slice(&body).unwrap())
    }
  }
}
