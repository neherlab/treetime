use app_backend::runs::errors::invalid;
use axum::body::Body;
use axum::http::{HeaderMap, Request};
use eyre::Report;

pub fn port_request(
  method: &str,
  url: &str,
  headers: impl IntoIterator<Item = (String, String)>,
  body: Option<Vec<u8>>,
) -> Result<Request<Body>, Report> {
  headers
    .into_iter()
    .fold(Request::builder().method(method).uri(url), |builder, (name, value)| {
      builder.header(name, value)
    })
    .body(body.map_or_else(Body::empty, Body::from))
    .map_err(|err| invalid(format!("the request `{method} {url}` is malformed: {err}")))
}

pub fn header_pairs(headers: &HeaderMap) -> Vec<(String, String)> {
  headers
    .iter()
    .map(|(name, value)| {
      (
        name.as_str().to_owned(),
        String::from_utf8_lossy(value.as_bytes()).into_owned(),
      )
    })
    .collect()
}
