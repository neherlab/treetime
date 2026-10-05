use app_commands::runs::errors::invalid;
use axum::body::Body;
use axum::http::{HeaderMap, Request};
use eyre::Report;
use napi::bindgen_prelude::Uint8Array;
use napi_derive::napi;

#[napi(object)]
pub struct PortRequest {
  pub seq: u32,
  pub method: String,
  pub url: String,
  pub headers: Vec<PortHeader>,
  pub body: Option<Uint8Array>,
}

impl PortRequest {
  pub fn into_http(self) -> Result<Request<Body>, Report> {
    let Self {
      method,
      url,
      headers,
      body,
      ..
    } = self;
    let builder = headers.into_iter().fold(
      Request::builder().method(method.as_str()).uri(url.as_str()),
      |builder, header| builder.header(header.name, header.value),
    );
    builder
      .body(body.map_or_else(Body::empty, |bytes| Body::from(bytes.to_vec())))
      .map_err(|err| invalid(format!("the request `{method} {url}` is malformed: {err}")))
  }
}

#[napi(object)]
pub struct PortHeader {
  pub name: String,
  pub value: String,
}

impl PortHeader {
  pub fn list(headers: &HeaderMap) -> Vec<Self> {
    headers
      .iter()
      .map(|(name, value)| Self {
        name: name.as_str().to_owned(),
        value: String::from_utf8_lossy(value.as_bytes()).into_owned(),
      })
      .collect()
  }
}

#[napi(discriminant = "kind", discriminant_case = "lowercase")]
pub enum PortMessage {
  Request { request: PortRequest },
  Abort { seq: u32 },
}

#[napi(discriminant = "kind", discriminant_case = "lowercase")]
pub enum PortReply {
  Head {
    seq: u32,
    status: u16,
    headers: Vec<PortHeader>,
  },
  Chunk {
    seq: u32,
    data: Uint8Array,
  },
  End {
    seq: u32,
  },
  Reset {
    seq: u32,
    message: String,
  },
}

impl PortReply {
  pub fn seq(&self) -> u32 {
    match self {
      Self::Head { seq, .. } | Self::Chunk { seq, .. } | Self::End { seq } | Self::Reset { seq, .. } => *seq,
    }
  }
}
