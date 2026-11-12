#[cfg(test)]
mod tests {
  use crate::__tests__::test_routes::tests::helpers::{app, body_json};
  use crate::error::panic_response;
  use axum::body::Body;
  use axum::http::Request;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown_path(      ("GET",    "/api/nope"),    (404, "not_found"))]
  #[case::unknown_root(      ("GET",    "/api"),         (404, "not_found"))]
  #[case::method_not_allowed(("DELETE", "/api/version"), (405, "method_not_allowed"))]
  #[trace]
  #[tokio::test]
  async fn test_errors_of_the_api_are_json(#[case] (method, path): (&str, &str), #[case] expected: (u16, &str)) {
    let test = app();
    let response = test
      .send(Request::builder().method(method).uri(path).body(Body::empty()).unwrap())
      .await;
    let (status, body) = body_json(response).await;
    assert_eq!((expected.0, json!(expected.1)), (status, body["code"].clone()));
  }

  #[tokio::test]
  async fn test_errors_panic_answers_the_json_internal_error() {
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
}
