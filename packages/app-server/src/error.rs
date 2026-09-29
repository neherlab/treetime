use app_commands::bridge::error::{ErrorCode, ErrorResponse};
use axum::Json;
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use eyre::Report;
use log::error;
use std::any::Any;

pub(crate) struct AppError(Report);

impl IntoResponse for AppError {
  fn into_response(self) -> Response {
    let body = ErrorResponse::from_report(&self.0);
    if body.code == ErrorCode::InternalError {
      error!("{:?}", self.0);
    }
    error_response(body)
  }
}

impl<E: Into<Report>> From<E> for AppError {
  fn from(err: E) -> Self {
    Self(err.into())
  }
}

pub(crate) fn error_response(body: ErrorResponse) -> Response {
  (http_status(body.code), Json(body)).into_response()
}

pub(crate) fn plain_error(code: ErrorCode, message: impl Into<String>) -> Response {
  error_response(ErrorResponse {
    code,
    message: message.into(),
    causes: vec![],
  })
}

pub(crate) fn panic_response(payload: &(dyn Any + Send)) -> Response {
  let body = ErrorResponse::from_panic(payload);
  error!("{}: {}", body.message, body.causes.join(": "));
  error_response(body)
}

fn http_status(code: ErrorCode) -> StatusCode {
  match code {
    ErrorCode::NotFound => StatusCode::NOT_FOUND,
    ErrorCode::UploadTooLarge => StatusCode::PAYLOAD_TOO_LARGE,
    ErrorCode::Conflict => StatusCode::CONFLICT,
    ErrorCode::InvalidRequest => StatusCode::BAD_REQUEST,
    ErrorCode::MethodNotAllowed => StatusCode::METHOD_NOT_ALLOWED,
    ErrorCode::Forbidden => StatusCode::FORBIDDEN,
    ErrorCode::Timeout => StatusCode::REQUEST_TIMEOUT,
    ErrorCode::InternalError => StatusCode::INTERNAL_SERVER_ERROR,
  }
}
