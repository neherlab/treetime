use crate::error::AppError;
use aide::OperationOutput;
use aide::generate::GenContext;
use aide::openapi::{MediaType, Operation, Response as ApiResponse, SchemaObject, StatusCode};
use app_commands::bridge::error::ErrorResponse;
use axum::http::{HeaderValue, header};
use axum::response::sse::{Event, KeepAlive, Sse};
use axum::response::{IntoResponse, Response};
use indexmap::IndexMap;
use schemars::{JsonSchema, Schema, json_schema};
use std::convert::Infallible;
use std::marker::PhantomData;
use std::pin::Pin;
use tokio_stream::{Stream, StreamExt as _};

const JSON: &str = "application/json";
const EVENT_STREAM: &str = "text/event-stream";
const OCTET_STREAM: &str = "application/octet-stream";
const ZIP: &str = "application/zip";

impl OperationOutput for AppError {
  type Inner = ErrorResponse;

  fn operation_response(ctx: &mut GenContext, _operation: &mut Operation) -> Option<ApiResponse> {
    let schema = ctx.schema.subschema_for::<ErrorResponse>();
    Some(content_response("The error, with its causes", JSON, schema))
  }

  fn inferred_responses(ctx: &mut GenContext, operation: &mut Operation) -> Vec<(Option<StatusCode>, ApiResponse)> {
    Self::operation_response(ctx, operation)
      .map(|response| (None, response))
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
    Sse::new(self.events.map(Ok::<_, Infallible>))
      .keep_alive(KeepAlive::default())
      .into_response()
  }
}

impl<T: JsonSchema> OperationOutput for TypedSse<T> {
  type Inner = T;

  fn operation_response(ctx: &mut GenContext, _operation: &mut Operation) -> Option<ApiResponse> {
    let schema = ctx.schema.subschema_for::<T>();
    Some(content_response(
      "Stream of server-sent events; the data of each event is one JSON item",
      EVENT_STREAM,
      schema,
    ))
  }

  fn inferred_responses(ctx: &mut GenContext, operation: &mut Operation) -> Vec<(Option<StatusCode>, ApiResponse)> {
    success(Self::operation_response(ctx, operation))
  }
}

pub(crate) struct FileContent {
  pub content_type: HeaderValue,
  pub bytes: Vec<u8>,
}

impl IntoResponse for FileContent {
  fn into_response(self) -> Response {
    ([(header::CONTENT_TYPE, self.content_type)], self.bytes).into_response()
  }
}

impl OperationOutput for FileContent {
  type Inner = Self;

  fn operation_response(_ctx: &mut GenContext, _operation: &mut Operation) -> Option<ApiResponse> {
    Some(content_response("Contents of the file", OCTET_STREAM, binary()))
  }

  fn inferred_responses(ctx: &mut GenContext, operation: &mut Operation) -> Vec<(Option<StatusCode>, ApiResponse)> {
    success(Self::operation_response(ctx, operation))
  }
}

pub(crate) struct ZipAttachment {
  pub disposition: HeaderValue,
  pub bytes: Vec<u8>,
}

impl IntoResponse for ZipAttachment {
  fn into_response(self) -> Response {
    (
      [
        (header::CONTENT_TYPE, HeaderValue::from_static(ZIP)),
        (header::CONTENT_DISPOSITION, self.disposition),
      ],
      self.bytes,
    )
      .into_response()
  }
}

impl OperationOutput for ZipAttachment {
  type Inner = Self;

  fn operation_response(_ctx: &mut GenContext, _operation: &mut Operation) -> Option<ApiResponse> {
    Some(content_response("Zip archive", ZIP, binary()))
  }

  fn inferred_responses(ctx: &mut GenContext, operation: &mut Operation) -> Vec<(Option<StatusCode>, ApiResponse)> {
    success(Self::operation_response(ctx, operation))
  }
}

fn success(response: Option<ApiResponse>) -> Vec<(Option<StatusCode>, ApiResponse)> {
  response
    .map(|response| (Some(StatusCode::Code(200)), response))
    .into_iter()
    .collect()
}

fn content_response(description: &str, content_type: &str, schema: Schema) -> ApiResponse {
  ApiResponse {
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
    ..ApiResponse::default()
  }
}

pub(crate) fn binary() -> Schema {
  json_schema!({ "type": "string", "format": "binary" })
}
