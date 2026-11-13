#[cfg(test)]
mod tests {
  use crate::port::{header_pairs, port_request};
  use axum::body::to_bytes;
  use axum::http::{HeaderMap, HeaderValue};
  use pretty_assertions::assert_eq;
  use treetime_utils::{assert_error, o};

  #[tokio::test]
  async fn test_port_request_carries_method_url_headers_and_body() {
    let body = vec![0xff, 0x00, 0xc3, 0x28];
    let request = port_request(
      "PUT",
      "/api/runs/abc/inputs/blob.bin?x=1",
      vec![
        (o!("content-type"), o!("application/octet-stream")),
        (o!("x-two"), o!("a")),
        (o!("x-two"), o!("b")),
      ],
      Some(body.clone()),
    )
    .unwrap();
    let (head, actual_body) = request.into_parts();
    let actual = (
      head.method.as_str().to_owned(),
      head.uri.to_string(),
      header_pairs(&head.headers),
      to_bytes(actual_body, usize::MAX).await.unwrap().to_vec(),
    );
    let expected = (
      o!("PUT"),
      o!("/api/runs/abc/inputs/blob.bin?x=1"),
      vec![
        (o!("content-type"), o!("application/octet-stream")),
        (o!("x-two"), o!("a")),
        (o!("x-two"), o!("b")),
      ],
      body,
    );
    assert_eq!(expected, actual);
  }

  #[tokio::test]
  async fn test_port_request_without_a_body_has_an_empty_body() {
    let request = port_request("GET", "/api/version", vec![], None).unwrap();
    let actual = to_bytes(request.into_body(), usize::MAX).await.unwrap().to_vec();
    assert_eq!(Vec::<u8>::new(), actual);
  }

  #[test]
  fn test_port_request_rejects_a_malformed_method() {
    assert_error!(
      port_request("NOT A METHOD", "/api/version", vec![], None),
      "the request `NOT A METHOD /api/version` is malformed: invalid HTTP method"
    );
  }

  #[test]
  fn test_port_request_rejects_a_malformed_header_name() {
    assert_error!(
      port_request("GET", "/api/version", vec![(o!("bad name"), o!("x"))], None),
      "the request `GET /api/version` is malformed: invalid HTTP header name"
    );
  }

  #[test]
  fn test_port_header_pairs_replaces_bytes_that_are_not_utf8() {
    let mut headers = HeaderMap::new();
    headers.append("x-text", HeaderValue::from_static("plain"));
    headers.append("x-bytes", HeaderValue::from_bytes(&[b'a', 0xff, b'b']).unwrap());
    let expected = vec![(o!("x-text"), o!("plain")), (o!("x-bytes"), o!("a\u{fffd}b"))];
    assert_eq!(expected, header_pairs(&headers));
  }
}
