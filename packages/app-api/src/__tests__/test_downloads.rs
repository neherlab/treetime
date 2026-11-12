#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, body_json, request, wait_for_status};
  use axum::body::Body;
  use axum::http::Request;
  use helpers::{archive_entries, finished_run, get, header_of};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use treetime_charts::convert::{svg_to_pdf, svg_to_png};
  use zip::CompressionMethod;

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
  async fn test_downloads_package_streams_stored_entries_of_the_run_folder_and_its_inputs() {
    let (test, id) = finished_run().await;
    let (_, files) = request(&test, "GET", &format!("/api/runs/{id}/files"), None).await;
    let expected = [
      "run.json",
      "run.yaml",
      "command.sh",
      "events.jsonl",
      "inputs/aln.fasta.xz",
      "inputs/metadata.tsv",
      "inputs/tree.nwk",
    ]
    .into_iter()
    .map(str::to_owned)
    .chain(
      files
        .as_array()
        .unwrap()
        .iter()
        .map(|file| format!("out/{}", file["path"].as_str().unwrap())),
    )
    .map(|name| (format!("treetime-{id}/{name}"), CompressionMethod::Stored))
    .sorted_by_key(|(name, _)| name.clone())
    .collect_vec();
    let response = get(&test, &format!("/api/runs/{id}/package"), &[]).await;
    let headers = (
      header_of(&response, "content-type"),
      header_of(&response, "content-disposition"),
    );
    let body = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
    assert_eq!(
      (
        (
          Some("application/zip".to_owned()),
          Some(format!(
            "attachment; filename=\"treetime-{id}.zip\"; filename*=UTF-8''treetime%2D{id}%2Ezip"
          ))
        ),
        expected
      ),
      (headers, archive_entries(body.to_vec()))
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_import_of_a_package_creates_a_run_with_the_same_results() {
    let (test, id) = finished_run().await;
    let response = get(&test, &format!("/api/runs/{id}/package"), &[]).await;
    let package = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();

    let response = test
      .send(
        Request::post("/api/runs/import")
          .header("content-type", "application/octet-stream")
          .body(Body::from(package))
          .unwrap(),
      )
      .await;
    let (status, imported) = body_json(response).await;
    let imported_id = imported["id"].as_str().unwrap();
    let (_, results) = request(&test, "GET", &format!("/api/runs/{id}/results"), None).await;
    let (_, imported_results) = request(&test, "GET", &format!("/api/runs/{imported_id}/results"), None).await;

    assert_eq!(
      (200, json!("ok"), true, results),
      (status, imported["status"].clone(), imported_id != id, imported_results)
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_import_of_a_file_that_is_not_a_zip_archive_is_an_invalid_request() {
    let test = app();
    let response = test
      .send(
        Request::post("/api/runs/import")
          .header("content-type", "application/octet-stream")
          .body(Body::from("not a zip archive"))
          .unwrap(),
      )
      .await;
    let (status, body) = body_json(response).await;
    assert_eq!((400, json!("invalid_request")), (status, body["code"].clone()));
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_run_file_charts_have_their_media_types() {
    let test = app();
    let (_, record) = test
      .run_with_data(
        "clock",
        json!({ "tree": "zika/20/tree.nwk", "metadata": "zika/20/metadata.tsv" }),
      )
      .await;
    let id = record["id"].as_str().unwrap().to_owned();
    wait_for_status(&test, &id, "ok").await;
    let out_dir = test.runs_dir.path().join(&id).join("out");
    let svg = fs::read_to_string(out_dir.join("clock.rtt.svg")).unwrap();
    fs::write(out_dir.join("clock.rtt.png"), svg_to_png(&svg).unwrap()).unwrap();
    fs::write(out_dir.join("clock.rtt.pdf"), svg_to_pdf(&svg).unwrap()).unwrap();
    let mut content_types = vec![];
    for path in ["clock.rtt.svg", "clock.rtt.png", "clock.rtt.pdf"] {
      let response = get(&test, &format!("/api/runs/{id}/file?path={path}"), &[]).await;
      content_types.push((response.status().as_u16(), header_of(&response, "content-type")));
    }
    assert_eq!(
      vec![
        (200, Some("image/svg+xml".to_owned())),
        (200, Some("image/png".to_owned())),
        (200, Some("application/pdf".to_owned())),
      ],
      content_types
    );
  }

  #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
  async fn test_downloads_package_of_an_unknown_run_is_not_found() {
    let test = app();
    let (status, body) = request(&test, "GET", "/api/runs/0123456789abcdef0123456789abcdef/package", None).await;
    assert_eq!((404, "not_found"), (status, body["code"].as_str().unwrap()));
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::{TestApp, app, timetree_config, wait_for_status};
    use axum::body::Body;
    use axum::http::Request;
    use axum::response::Response;
    use itertools::Itertools;
    use std::io::Cursor;
    use zip::{CompressionMethod, ZipArchive};

    pub(super) async fn finished_run() -> (TestApp, String) {
      let test = app();
      let (_, record) = test.run_with_data("timetree", timetree_config()).await;
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
          (file.name().unwrap().into_owned(), file.compression())
        })
        .sorted_by_key(|(name, _)| name.clone())
        .collect_vec()
    }
  }
}
