use crate::api::response::binary;
use crate::error::AppError;
use aide::generate::GenContext;
use aide::openapi::{
  MediaType, Operation, RequestBody, Response as ApiResponse, SchemaObject, StatusCode as ApiStatusCode,
};
use aide::operation::set_body;
use aide::{OperationInput, OperationOutput};
use app_commands::runs::errors::invalid;
use axum::Json;
use axum::body::{Body, Bytes};
use axum::extract::{FromRequest, FromRequestParts, Path, Query, RawPathParams, Request};
use axum::http::header::CONTENT_TYPE;
use axum::http::request::Parts;
use axum::http::{HeaderMap, HeaderValue};
use axum::response::{IntoResponse, Response};
use deser::Serialize;
use deser::de::DeserializeOwned;
use deser_path::PathLayer;
use deser_value::{Kind, Map, Value, from_value};
use indexmap::IndexMap;
use mime::Mime;
use schemars::JsonSchema;
use std::convert::Infallible;
use std::future::{Future, ready};
use treetime_utils::io::json::{JsonPretty, json_write};

pub(crate) struct ApiJson<T>(pub T);

impl<T, S> FromRequest<S> for ApiJson<T>
where
  T: DeserializeOwned,
  S: Send + Sync,
{
  type Rejection = AppError;

  async fn from_request(request: Request, state: &S) -> Result<Self, AppError> {
    if !is_json(request.headers()) {
      return Err(invalid("Expected request with `Content-Type: application/json`").into());
    }
    let body = Bytes::from_request(request, state)
      .await
      .map_err(|rejection| invalid(rejection.body_text()))?;
    deser_json::Deserializer::from_slice(&body)
      .deserialize_with(|driver| driver.push_layer(PathLayer::new()))
      .map(Self)
      .map_err(|err| invalid(format!("Failed to parse the request body as JSON: {err}")).into())
  }
}

impl<T: JsonSchema> OperationInput for ApiJson<T> {
  fn operation_input(ctx: &mut GenContext, operation: &mut Operation) {
    Json::<T>::operation_input(ctx, operation);
  }
}

impl<T: Serialize> IntoResponse for ApiJson<T> {
  fn into_response(self) -> Response {
    let mut body = Vec::new();
    match json_write(&mut body, &self.0, JsonPretty(false)) {
      Ok(()) => (
        [(CONTENT_TYPE, HeaderValue::from_static(mime::APPLICATION_JSON.as_ref()))],
        body,
      )
        .into_response(),
      Err(err) => AppError::from(err).into_response(),
    }
  }
}

impl<T: JsonSchema> OperationOutput for ApiJson<T> {
  type Inner = T;

  fn operation_response(ctx: &mut GenContext, operation: &mut Operation) -> Option<ApiResponse> {
    Json::<T>::operation_response(ctx, operation)
  }

  fn inferred_responses(ctx: &mut GenContext, operation: &mut Operation) -> Vec<(Option<ApiStatusCode>, ApiResponse)> {
    Json::<T>::inferred_responses(ctx, operation)
  }
}

pub(crate) struct ApiPath<T>(pub T);

impl<T, S> FromRequestParts<S> for ApiPath<T>
where
  T: DeserializeOwned + Send,
  S: Send + Sync,
{
  type Rejection = AppError;

  async fn from_request_parts(parts: &mut Parts, state: &S) -> Result<Self, AppError> {
    let params = RawPathParams::from_request_parts(parts, state)
      .await
      .map_err(|rejection| invalid(rejection.body_text()))?;
    let params: Map = params
      .iter()
      .map(|(key, value)| (key, Value::from(Kind::Lexical(value.to_owned()))))
      .collect();
    from_value(&Value::from(params))
      .map(Self)
      .map_err(|err| invalid(format!("Invalid URL: {err}")).into())
  }
}

impl<T: JsonSchema> OperationInput for ApiPath<T> {
  fn operation_input(ctx: &mut GenContext, operation: &mut Operation) {
    Path::<T>::operation_input(ctx, operation);
  }
}

pub(crate) struct ApiQuery<T>(pub T);

impl<T, S> FromRequestParts<S> for ApiQuery<T>
where
  T: DeserializeOwned,
  S: Send + Sync,
{
  type Rejection = AppError;

  fn from_request_parts(parts: &mut Parts, _state: &S) -> impl Future<Output = Result<Self, AppError>> {
    ready(
      deser_urlencoded::Deserializer::from_str(parts.uri.query().unwrap_or_default())
        .deserialize_with(|driver| driver.push_layer(PathLayer::new()))
        .map(Self)
        .map_err(|err| invalid(format!("Failed to deserialize query string: {err}")).into()),
    )
  }
}

impl<T: JsonSchema> OperationInput for ApiQuery<T> {
  fn operation_input(ctx: &mut GenContext, operation: &mut Operation) {
    Query::<T>::operation_input(ctx, operation);
  }
}

pub(crate) struct OctetStream(pub Body);

impl<S: Send + Sync> FromRequest<S> for OctetStream {
  type Rejection = Infallible;

  fn from_request(request: Request, _state: &S) -> impl Future<Output = Result<Self, Infallible>> {
    ready(Ok(Self(request.into_body())))
  }
}

impl OperationInput for OctetStream {
  fn operation_input(ctx: &mut GenContext, operation: &mut Operation) {
    set_body(
      ctx,
      operation,
      RequestBody {
        description: None,
        content: IndexMap::from_iter([(
          "application/octet-stream".to_owned(),
          MediaType {
            schema: Some(SchemaObject {
              json_schema: binary(),
              example: None,
              external_docs: None,
            }),
            ..MediaType::default()
          },
        )]),
        required: true,
        extensions: IndexMap::default(),
      },
    );
  }
}

fn is_json(headers: &HeaderMap) -> bool {
  headers
    .get(CONTENT_TYPE)
    .and_then(|value| value.to_str().ok())
    .and_then(|value| value.parse::<Mime>().ok())
    .is_some_and(|mime| {
      mime.type_() == mime::APPLICATION
        && (mime.subtype() == mime::JSON || mime.suffix().is_some_and(|suffix| suffix == mime::JSON))
    })
}
