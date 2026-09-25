use app_commands::bridge::error::{ErrorCode, ErrorResponse};
use axum::Json;
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use eyre::Report;
use log::error;

pub(crate) struct AppError(Report);

impl IntoResponse for AppError {
  fn into_response(self) -> Response {
    let body = ErrorResponse::from_report(&self.0);
    let status = http_status(body.code);
    if status == StatusCode::INTERNAL_SERVER_ERROR {
      error!("{:?}", self.0);
    }
    (status, Json(body)).into_response()
  }
}

impl<E: Into<Report>> From<E> for AppError {
  fn from(err: E) -> Self {
    Self(err.into())
  }
}

fn http_status(code: ErrorCode) -> StatusCode {
  match code {
    ErrorCode::NotFound => StatusCode::NOT_FOUND,
    ErrorCode::UploadTooLarge => StatusCode::PAYLOAD_TOO_LARGE,
    ErrorCode::Conflict => StatusCode::CONFLICT,
    ErrorCode::InvalidRequest => StatusCode::BAD_REQUEST,
    ErrorCode::InternalError => StatusCode::INTERNAL_SERVER_ERROR,
  }
}
