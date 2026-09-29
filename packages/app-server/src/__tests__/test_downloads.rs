#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, create_deferred, request};
  use helpers::{archive_entries, finished_run, get, header_of};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use zip::CompressionMethod;

  #[rustfmt::skip]
  #[rstest]
  #[case::results(    "/api/runs/ID/results")]
  #[case::auspice(    "/api/runs/ID/auspice")]
  #[case::comparison( "/api/runs/ID/compare/ID")]
  #[trace]
  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_finished_results_revalidate_with_their_etag(#[case] template: &str) {
    let (test, id) = finished_run().await;
    let path = template.replace("ID", &id);
    let first = get(&test, &path, &[]).await;
    let etag = header_of(&first, "etag").unwrap();
    let second = get(&test, &path, &[("if-none-match", &etag)]).await;
    let other = get(&test, &path, &[("if-none-match", "W/\"0000000000000000\"")]).await;
    let body = axum::body::to_bytes(second.into_body(), usize::MAX).await.unwrap();
    assert_eq!(
      (
        (200, true, Some("private, no-cache".to_owned())),
        (304, true),
        200
      ),
      (
        (
          first.status().as_u16(),
          etag.starts_with("W/\""),
          header_of(&first, "cache-control")
        ),
        (304, body.is_empty()),
        other.status().as_u16()
      )
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_results_of_an_unfinished_run_have_no_etag() {
    let (test, _) = finished_run().await;
    let id = create_deferred(&test).await;
    let response = get(&test, &format!("/api/runs/{id}/results"), &[]).await;
    assert_eq!(
      (None, Some("no-store".to_owned())),
      (header_of(&response, "etag"), header_of(&response, "cache-control"))
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_run_file_serves_byte_ranges_as_an_attachment() {
    let (test, id) = finished_run().await;
    let response = get(
      &test,
      &format!("/api/runs/{id}/file?path=timetree.auspice.json"),
      &[("range", "bytes=0-3")],
    )
    .await;
    let status = response.status().as_u16();
    let headers = (
      header_of(&response, "content-type"),
      header_of(&response, "content-disposition"),
      header_of(&response, "content-range").map(|range| range.starts_with("bytes 0-3/")),
    );
    let body = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
    assert_eq!(
      (
        206,
        (
          Some("application/json".to_owned()),
          Some("attachment; filename=\"timetree.auspice.json\"; filename*=UTF-8''timetree%2Eauspice%2Ejson".to_owned()),
          Some(true)
        ),
        4
      ),
      (status, headers, body.len())
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_archive_streams_stored_entries_of_every_output_file() {
    let (test, id) = finished_run().await;
    let (_, files) = request(&test, "GET", &format!("/api/runs/{id}/files"), None).await;
    let expected = files
      .as_array()
      .unwrap()
      .iter()
      .map(|file| {
        (
          format!("{id}/{}", file["path"].as_str().unwrap()),
          CompressionMethod::Stored,
        )
      })
      .sorted_by_key(|(name, _)| name.clone())
      .collect_vec();
    let response = get(&test, &format!("/api/runs/{id}/archive"), &[]).await;
    let content_type = header_of(&response, "content-type");
    let body = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
    assert_eq!(
      (Some("application/zip".to_owned()), expected),
      (content_type, archive_entries(body.to_vec()))
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_archive_of_an_unknown_run_is_not_found() {
    let test = app();
    let (status, body) = request(&test, "GET", "/api/runs/0123456789abcdef0123456789abcdef/archive", None).await;
    assert_eq!((404, "not_found"), (status, body["code"].as_str().unwrap()));
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::{TestApp, app, request, timetree_config, wait_for_status};
    use axum::body::Body;
    use axum::http::Request;
    use axum::response::Response;
    use itertools::Itertools;
    use serde_json::json;
    use std::io::Cursor;
    use zip::{CompressionMethod, ZipArchive};

    pub(super) async fn finished_run() -> (TestApp, String) {
      let test = app();
      let (_, record) = request(
        &test,
        "POST",
        "/api/runs",
        Some(json!({ "command": "timetree", "config": timetree_config() })),
      )
      .await;
      let id = record["id"].as_str().unwrap().to_owned();
      wait_for_status(&test, &id, "ok").await;
      (test, id)
    }

    pub(super) async fn get(test: &TestApp, path: &str, headers: &[(&str, &str)]) -> Response {
      let request = headers.iter().fold(Request::get(path), |request, (name, value)| {
        request.header(*name, *value)
      });
      test.send(request.body(Body::empty()).unwrap()).await
    }

    pub(super) fn header_of(response: &Response, name: &str) -> Option<String> {
      response
        .headers()
        .get(name)
        .map(|value| value.to_str().unwrap().to_owned())
    }

    pub(super) fn archive_entries(bytes: Vec<u8>) -> Vec<(String, CompressionMethod)> {
      let mut archive = ZipArchive::new(Cursor::new(bytes)).unwrap();
      (0..archive.len())
        .map(|index| {
          let file = archive.by_index(index).unwrap();
          (file.name().to_owned(), file.compression())
        })
        .sorted_by_key(|(name, _)| name.clone())
        .collect_vec()
    }
  }
}
