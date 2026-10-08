use app_server::port::{header_pairs, port_request};
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
    port_request(
      &method,
      &url,
      headers.into_iter().map(|header| (header.name, header.value)),
      body.map(|bytes| bytes.to_vec()),
    )
  }
}

#[napi(object)]
pub struct PortHeader {
  pub name: String,
  pub value: String,
}

impl PortHeader {
  pub fn list(headers: &HeaderMap) -> Vec<Self> {
    header_pairs(headers)
      .into_iter()
      .map(|(name, value)| Self { name, value })
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
