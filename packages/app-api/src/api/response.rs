use crate::api::extract::ApiJson;
use crate::error::AppError;
use aide::OperationOutput;
use aide::generate::GenContext;
use aide::openapi::{self, MediaType, Operation, SchemaObject};
use app_backend::bridge::error::ErrorResponse;
use axum::body::Body;
use axum::http::{HeaderValue, StatusCode, header};
use axum::response::sse::{Event, Sse};
use axum::response::{IntoResponse, Response};
use deser::Serialize;
use indexmap::IndexMap;
use schemars::{JsonSchema, Schema, json_schema};
use std::convert::Infallible;
use std::marker::PhantomData;
use std::pin::Pin;
use tokio_stream::{Stream, StreamExt as _};

#[cfg(not(target_family = "wasm"))]
use axum::response::sse::KeepAlive;

const JSON: &str = "application/json";
const EVENT_STREAM: &str = "text/event-stream";
const OCTET_STREAM: &str = "application/octet-stream";
const ZIP: &str = "application/zip";

impl OperationOutput for AppError {
  type Inner = ErrorResponse;

  fn operation_response(ctx: &mut GenContext, _operation: &mut Operation) -> Option<openapi::Response> {
    let schema = ctx.schema.subschema_for::<ErrorResponse>();
    Some(content_response("The error, with its causes", JSON, schema))
  }

  fn inferred_responses(
    ctx: &mut GenContext,
    operation: &mut Operation,
  ) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
    Self::operation_response(ctx, operation)
      .map(|response| (None, response))
      .into_iter()
      .collect()
  }
}

pub(crate) struct Accepted<T>(pub T);

impl<T: Serialize> IntoResponse for Accepted<T> {
  fn into_response(self) -> Response {
    (StatusCode::ACCEPTED, ApiJson(self.0)).into_response()
  }
}

impl<T: JsonSchema> OperationOutput for Accepted<T> {
  type Inner = T;

  fn operation_response(ctx: &mut GenContext, _operation: &mut Operation) -> Option<openapi::Response> {
    let schema = ctx.schema.subschema_for::<T>();
    Some(content_response(
      "The request was accepted and runs in the background",
      JSON,
      schema,
    ))
  }

  fn inferred_responses(
    ctx: &mut GenContext,
    operation: &mut Operation,
  ) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
    Self::operation_response(ctx, operation)
      .map(|response| (Some(openapi::StatusCode::Code(StatusCode::ACCEPTED.as_u16())), response))
      .into_iter()
      .collect()
  }
}

pub(crate) struct TypedSse<T> {
  events: Pin<Box<dyn Stream<Item = Event> + Send>>,
  item: PhantomData<fn() -> T>,
}

impl<T: 'static> TypedSse<T> {
  pub(crate) fn new(items: impl Stream<Item = T> + Send + 'static, to_event: fn(&T) -> Event) -> Self {
    Self {
      events: Box::pin(items.map(move |item| to_event(&item))),
      item: PhantomData,
    }
  }
}

impl<T> IntoResponse for TypedSse<T> {
  fn into_response(self) -> Response {
    let sse = Sse::new(self.events.map(Ok::<_, Infallible>));
    #[cfg(not(target_family = "wasm"))]
    let sse = sse.keep_alive(KeepAlive::default());
    sse.into_response()
  }
}

impl<T: JsonSchema> OperationOutput for TypedSse<T> {
  type Inner = T;

  fn operation_response(ctx: &mut GenContext, _operation: &mut Operation) -> Option<openapi::Response> {
    let schema = ctx.schema.subschema_for::<T>();
    Some(content_response(
      "Stream of server-sent events; the data of each event is one JSON item",
      EVENT_STREAM,
      schema,
    ))
  }

  fn inferred_responses(
    ctx: &mut GenContext,
    operation: &mut Operation,
  ) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
    success(Self::operation_response(ctx, operation))
  }
}

pub(crate) struct FileContent(pub Response);

impl IntoResponse for FileContent {
  fn into_response(self) -> Response {
    self.0
  }
}

impl OperationOutput for FileContent {
  type Inner = Self;

  fn operation_response(_ctx: &mut GenContext, _operation: &mut Operation) -> Option<openapi::Response> {
    Some(content_response("Contents of the file", OCTET_STREAM, binary()))
  }

  fn inferred_responses(
    ctx: &mut GenContext,
    operation: &mut Operation,
  ) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
    success(Self::operation_response(ctx, operation))
  }
}

pub(crate) struct ZipAttachment {
  pub disposition: HeaderValue,
  pub body: Body,
}

impl IntoResponse for ZipAttachment {
  fn into_response(self) -> Response {
    (
      [
        (header::CONTENT_TYPE, HeaderValue::from_static(ZIP)),
        (header::CONTENT_DISPOSITION, self.disposition),
      ],
      self.body,
    )
      .into_response()
  }
}

impl OperationOutput for ZipAttachment {
  type Inner = Self;

  fn operation_response(_ctx: &mut GenContext, _operation: &mut Operation) -> Option<openapi::Response> {
    Some(content_response("Zip archive", ZIP, binary()))
  }

  fn inferred_responses(
    ctx: &mut GenContext,
    operation: &mut Operation,
  ) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
    success(Self::operation_response(ctx, operation))
  }
}

fn success(response: Option<openapi::Response>) -> Vec<(Option<openapi::StatusCode>, openapi::Response)> {
  response
    .map(|response| (Some(openapi::StatusCode::Code(StatusCode::OK.as_u16())), response))
    .into_iter()
    .collect()
}

fn content_response(description: &str, content_type: &str, schema: Schema) -> openapi::Response {
  openapi::Response {
    description: description.to_owned(),
    content: IndexMap::from_iter([(
      content_type.to_owned(),
      MediaType {
        schema: Some(SchemaObject {
          json_schema: schema,
          example: None,
          external_docs: None,
        }),
        ..MediaType::default()
      },
    )]),
    ..openapi::Response::default()
  }
}

pub(crate) fn binary() -> Schema {
  json_schema!({ "type": "string", "format": "binary" })
}
