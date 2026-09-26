use crate::api::response::binary;
use crate::error::AppError;
use aide::OperationInput;
use aide::generate::GenContext;
use aide::openapi::{MediaType, Operation, RequestBody, SchemaObject};
use aide::operation::set_body;
use app_commands::runs::errors::invalid;
use axum::Json;
use axum::body::Body;
use axum::extract::{FromRequest, FromRequestParts, Path, Query, Request};
use axum::http::request::Parts;
use indexmap::IndexMap;
use schemars::JsonSchema;
use serde::de::DeserializeOwned;
use std::convert::Infallible;
use std::future::{Future, ready};

pub(crate) struct ApiJson<T>(pub T);

impl<T, S> FromRequest<S> for ApiJson<T>
where
  T: DeserializeOwned,
  S: Send + Sync,
{
  type Rejection = AppError;

  async fn from_request(request: Request, state: &S) -> Result<Self, AppError> {
    let Json(value) = Json::<T>::from_request(request, state)
      .await
      .map_err(|rejection| invalid(rejection.body_text()))?;
    Ok(Self(value))
  }
}

impl<T: JsonSchema> OperationInput for ApiJson<T> {
  fn operation_input(ctx: &mut GenContext, operation: &mut Operation) {
    Json::<T>::operation_input(ctx, operation);
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
    let Path(value) = Path::<T>::from_request_parts(parts, state)
      .await
      .map_err(|rejection| invalid(rejection.body_text()))?;
    Ok(Self(value))
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

  async fn from_request_parts(parts: &mut Parts, state: &S) -> Result<Self, AppError> {
    let Query(value) = Query::<T>::from_request_parts(parts, state)
      .await
      .map_err(|rejection| invalid(rejection.body_text()))?;
    Ok(Self(value))
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
