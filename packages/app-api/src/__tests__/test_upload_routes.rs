#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app_with_upload_limit, body_json, create_deferred, request};
  use helpers::{chunked_upload, stored_files};
  use percent_encoding::{NON_ALPHANUMERIC, utf8_percent_encode};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use std::path::PathBuf;
  use treetime_utils::o;

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_upload_routes_store_a_body_that_arrives_in_chunks() {
    let test = app_with_upload_limit(1000);
    let id = create_deferred(&test).await;

    let (status, uploaded) = body_json(
      test
        .send(chunked_upload(&id, "aln.fasta", &[b">A\n", b"AC", b"GT\n"]))
        .await,
    )
    .await;

    assert_eq!(
      (200, json!(8), o!(">A\nACGT\n")),
      (
        status,
        uploaded["size"].clone(),
        fs::read_to_string(uploaded["path"].as_str().unwrap()).unwrap()
      )
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_upload_routes_refuse_a_stream_over_the_limit_and_keep_no_file() {
    let test = app_with_upload_limit(1000);
    let id = create_deferred(&test).await;

    let (status, body) = body_json(
      test
        .send(chunked_upload(&id, "a.nwk", &[&[b'x'; 600], &[b'x'; 600]]))
        .await,
    )
    .await;

    assert_eq!(
      (
        413,
        json!("the inputs of a run are limited to 1000 bytes in total; `a.nwk` does not fit"),
        Vec::<PathBuf>::new(),
      ),
      (
        status,
        body["message"].clone(),
        stored_files(&test.runs_dir.path().join(".inputs"))
      ),
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_upload_routes_report_the_allowance_without_the_replaced_uploads() {
    let test = app_with_upload_limit(1000);
    let id = create_deferred(&test).await;
    let (_, first) = test.upload(&id, "tree.nwk", vec![b'x'; 300]).await;
    let first_path = utf8_percent_encode(first["path"].as_str().unwrap(), NON_ALPHANUMERIC).to_string();

    let all = request(&test, "GET", &format!("/api/runs/{id}/inputs"), None).await;
    let replacing = request(
      &test,
      "GET",
      &format!("/api/runs/{id}/inputs?replaces={first_path}"),
      None,
    )
    .await;

    assert_eq!(
      (
        (
          200,
          json!({ "remaining": 700, "limit": 1000, "limit_text": "1000 bytes" })
        ),
        (
          200,
          json!({ "remaining": 1000, "limit": 1000, "limit_text": "1000 bytes" })
        ),
      ),
      (all, replacing)
    );
  }

  mod helpers {
    use axum::body::{Body, Bytes};
    use axum::http::Request;
    use itertools::Itertools;
    use std::convert::Infallible;
    use std::fs;
    use std::path::{Path, PathBuf};
    use tokio_stream::iter;

    pub(super) fn chunked_upload(id: &str, name: &str, chunks: &[&[u8]]) -> Request<Body> {
      let chunks = chunks
        .iter()
        .map(|chunk| Ok::<_, Infallible>(Bytes::copy_from_slice(chunk)))
        .collect_vec();
      Request::put(format!("/api/runs/{id}/inputs/{name}"))
        .body(Body::from_stream(iter(chunks)))
        .unwrap()
    }

    pub(super) fn stored_files(dir: &Path) -> Vec<PathBuf> {
      fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .collect_vec()
    }
  }
}
