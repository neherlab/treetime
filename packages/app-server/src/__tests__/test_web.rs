#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::body_json;
  use crate::error::panic_response;
  use axum::body::Body;
  use axum::http::Request;
  use helpers::{header_of, headers_of, web_app};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use std::fs;
  use std::time::Duration;

  const IMMUTABLE: &str = "public, max-age=31536000, immutable";

  #[rustfmt::skip]
  #[rstest]
  #[case::asset_brotli(       ("/assets/app-0123abcd.js", "*/*",       "br, gzip"), (200, Some("br"),   Some(IMMUTABLE)))]
  #[case::asset_gzip(         ("/assets/app-0123abcd.js", "*/*",       "gzip"),     (200, Some("gzip"), Some(IMMUTABLE)))]
  #[case::asset_identity(     ("/assets/app-0123abcd.js", "*/*",       "identity"), (200, None,         Some(IMMUTABLE)))]
  #[case::asset_missing(      ("/assets/app-deadbeef.js", "*/*",       "gzip"),     (404, None,         None))]
  #[case::index(              ("/",                       "text/html", "identity"), (200, None,         Some("no-cache")))]
  #[case::page_navigation(    ("/runs/abc",               "text/html", "identity"), (200, None,         Some("no-cache")))]
  #[case::page_navigation_br( ("/runs/abc",               "text/html", "br"),       (200, Some("br"),   Some("no-cache")))]
  #[case::not_a_navigation(   ("/runs/abc.js",            "*/*",       "identity"), (404, None,         None))]
  #[case::public_file(        ("/favicon.svg",            "image/*",   "identity"), (200, None,         Some("no-cache")))]
  #[case::unknown_api_path(   ("/api/nope",               "text/html", "identity"), (404, None,         Some("no-store")))]
  #[trace]
  #[tokio::test]
  async fn test_web_static_files_are_cached_by_kind(
    #[case] (path, accept, encoding): (&str, &str, &str),
    #[case] expected: (u16, Option<&str>, Option<&str>),
  ) {
    let app = web_app();
    let response = app.get(path, &[("accept", accept), ("accept-encoding", encoding)]).await;
    let actual = (
      response.status().as_u16(),
      header_of(&response, "content-encoding"),
      header_of(&response, "cache-control"),
    );
    assert_eq!(expected, (actual.0, actual.1.as_deref(), actual.2.as_deref()));
  }

  #[tokio::test]
  async fn test_web_navigation_fallback_answers_the_index_page() {
    let app = web_app();
    let response = app.get("/runs/abc", &[("accept", "text/html")]).await;
    let body = axum::body::to_bytes(response.into_body(), usize::MAX).await.unwrap();
    let index = fs::read(app.dist.path().join("index.html")).unwrap();
    assert_eq!(index, body.to_vec());
  }

  #[tokio::test]
  async fn test_web_precompressed_assets_vary_by_accept_encoding() {
    let app = web_app();
    let response = app.get("/assets/app-0123abcd.js", &[("accept-encoding", "br")]).await;
    assert_eq!(Some("accept-encoding".to_owned()), header_of(&response, "vary"));
  }

  #[tokio::test]
  async fn test_web_assets_revalidate_with_their_etag() {
    let app = web_app();
    let first = app.get("/assets/app-0123abcd.js", &[]).await;
    let etag = header_of(&first, "etag").unwrap();
    let second = app.get("/assets/app-0123abcd.js", &[("if-none-match", &etag)]).await;
    assert_eq!(304, second.status().as_u16());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown_path(      ("GET",    "/api/nope"),   (404, "not_found"))]
  #[case::unknown_root(      ("GET",    "/api"),        (404, "not_found"))]
  #[case::method_not_allowed(("DELETE", "/api/health"), (405, "method_not_allowed"))]
  #[trace]
  #[tokio::test]
  async fn test_web_api_errors_are_json(#[case] (method, path): (&str, &str), #[case] expected: (u16, &str)) {
    let app = web_app();
    let response = app
      .test
      .send(Request::builder().method(method).uri(path).body(Body::empty()).unwrap())
      .await;
    let (status, body) = body_json(response).await;
    assert_eq!((expected.0, json!(expected.1)), (status, body["code"].clone()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::localhost_with_port(  Some("localhost:5173"),           200)]
  #[case::loopback_v4(          Some("127.0.0.1:3100"),           200)]
  #[case::loopback_v6(          Some("[::1]:3100"),               200)]
  #[case::localhost_subdomain(  Some("feat.treetime.localhost"),  200)]
  #[case::configured_host(      Some("treetime.example"),         200)]
  #[case::configured_host_case( Some("TreeTime.Example:443"),     200)]
  #[case::other_host(           Some("evil.example"),             403)]
  #[case::rebinding_suffix(     Some("treetime.example.evil"),    403)]
  #[case::no_host(              None,                             200)]
  #[trace]
  #[tokio::test]
  async fn test_web_hosts_outside_the_allowlist_are_forbidden(#[case] host: Option<&str>, #[case] expected: u16) {
    let app = web_app();
    let headers = host.map(|host| vec![("host", host)]).unwrap_or_default();
    let response = app.get("/api/health", &headers).await;
    assert_eq!(expected, response.status().as_u16());
  }

  #[tokio::test]
  async fn test_web_forbidden_host_answers_the_json_error() {
    let app = web_app();
    let response = app.get("/api/health", &[("host", "evil.example")]).await;
    let (status, body) = body_json(response).await;
    assert_eq!(
      (
        403,
        json!({
          "code": "forbidden",
          "message": "the server does not answer requests for host `evil.example`",
          "causes": [],
        })
      ),
      (status, body)
    );
  }

  #[tokio::test]
  async fn test_web_responses_carry_security_headers() {
    let app = web_app();
    let response = app.get("/api/health", &[]).await;
    assert_eq!(
      vec![
        Some("no-store".to_owned()),
        Some("nosniff".to_owned()),
        Some("strict-origin-when-cross-origin".to_owned()),
        Some("same-origin".to_owned()),
        Some("same-origin".to_owned()),
        Some("frame-ancestors 'none'".to_owned()),
        None,
      ],
      headers_of(
        &response,
        &[
          "cache-control",
          "x-content-type-options",
          "referrer-policy",
          "cross-origin-opener-policy",
          "cross-origin-resource-policy",
          "content-security-policy",
          "access-control-allow-origin",
        ]
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zstd_preferred("gzip, zstd", Some("zstd"))]
  #[case::gzip_only(     "gzip",       Some("gzip"))]
  #[case::brotli_only(   "br",         None)]
  #[case::identity(      "identity",   None)]
  #[trace]
  #[tokio::test]
  async fn test_web_api_responses_are_compressed(#[case] encoding: &str, #[case] expected: Option<&str>) {
    let app = web_app();
    let response = app.get("/api/openapi.json", &[("accept-encoding", encoding)]).await;
    assert_eq!(expected, header_of(&response, "content-encoding").as_deref());
  }

  #[tokio::test]
  async fn test_web_event_streams_are_not_compressed() {
    let app = web_app();
    let response = app.get("/api/events", &[("accept-encoding", "zstd, gzip")]).await;
    assert_eq!(
      vec![Some("text/event-stream".to_owned()), None],
      headers_of(&response, &["content-type", "content-encoding"])
    );
  }

  #[tokio::test]
  async fn test_web_event_streams_end_on_shutdown() {
    let app = web_app();
    let response = app.get("/api/events", &[]).await;
    app.shutdown.cancel();
    let body = tokio::time::timeout(
      Duration::from_secs(10),
      axum::body::to_bytes(response.into_body(), usize::MAX),
    )
    .await;
    assert!(matches!(body, Ok(Ok(_))), "the stream ends after shutdown");
  }

  #[tokio::test]
  async fn test_web_panic_answers_the_json_internal_error() {
    let (status, body) = body_json(panic_response(&"index out of bounds")).await;
    assert_eq!(
      (
        500,
        json!({
          "code": "internal_error",
          "message": "the back end stopped the operation after an internal error",
          "causes": ["index out of bounds"],
        })
      ),
      (status, body)
    );
  }

  mod helpers {
    use crate::__tests__::test_routes::tests::helpers::{TestApp, app_with};
    use crate::state::DEFAULT_MAX_UPLOAD_SIZE;
    use crate::web::WebOptions;
    use axum::body::Body;
    use axum::http::Request;
    use axum::response::Response;
    use std::fs;
    use tempfile::{TempDir, tempdir};
    use tokio_util::sync::CancellationToken;

    pub(super) struct WebApp {
      pub test: TestApp,
      pub shutdown: CancellationToken,
      pub dist: TempDir,
    }

    impl WebApp {
      pub(super) async fn get(&self, path: &str, headers: &[(&str, &str)]) -> Response {
        let request = headers.iter().fold(Request::get(path), |request, (name, value)| {
          request.header(*name, *value)
        });
        self.test.send(request.body(Body::empty()).unwrap()).await
      }
    }

    pub(super) fn web_app() -> WebApp {
      let dist = tempdir().unwrap();
      let assets = dist.path().join("assets");
      fs::create_dir_all(&assets).unwrap();
      fs::write(assets.join("app-0123abcd.js"), "console.log('app');".repeat(20)).unwrap();
      fs::write(assets.join("app-0123abcd.js.gz"), "gzip bytes").unwrap();
      fs::write(assets.join("app-0123abcd.js.br"), "brotli bytes").unwrap();
      fs::write(dist.path().join("index.html"), "<!doctype html><title>TreeTime</title>").unwrap();
      fs::write(dist.path().join("index.html.br"), "brotli index").unwrap();
      fs::write(dist.path().join("favicon.svg"), "<svg/>").unwrap();
      let shutdown = CancellationToken::new();
      let options = WebOptions {
        static_dir: Some(dist.path().to_path_buf()),
        allowed_hosts: vec!["treetime.example".to_owned()],
      };
      WebApp {
        test: app_with(DEFAULT_MAX_UPLOAD_SIZE, &shutdown, &options),
        shutdown,
        dist,
      }
    }

    pub(super) fn header_of(response: &Response, name: &str) -> Option<String> {
      response
        .headers()
        .get(name)
        .map(|value| value.to_str().unwrap().to_owned())
    }

    pub(super) fn headers_of(response: &Response, names: &[&str]) -> Vec<Option<String>> {
      names.iter().map(|name| header_of(response, name)).collect()
    }
  }
}
