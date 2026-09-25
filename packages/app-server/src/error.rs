use app_commands::runs::errors::{InvalidRunRequest, RunConflict, RunNotFound, UploadTooLarge};
use axum::Json;
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use eyre::Report;
use log::error;
use serde_json::json;

pub(crate) struct AppError(Report);

impl AppError {
  fn status(&self) -> (StatusCode, &'static str) {
    let report = &self.0;
    if report.downcast_ref::<RunNotFound>().is_some() {
      (StatusCode::NOT_FOUND, "not_found")
    } else if report.downcast_ref::<UploadTooLarge>().is_some() {
      (StatusCode::PAYLOAD_TOO_LARGE, "upload_too_large")
    } else if report.downcast_ref::<RunConflict>().is_some() {
      (StatusCode::CONFLICT, "conflict")
    } else if report.downcast_ref::<InvalidRunRequest>().is_some()
      || report.downcast_ref::<serde_json::Error>().is_some()
    {
      (StatusCode::BAD_REQUEST, "invalid_request")
    } else {
      (StatusCode::INTERNAL_SERVER_ERROR, "internal_error")
    }
  }
}

impl IntoResponse for AppError {
  fn into_response(self) -> Response {
    let (status, code) = self.status();
    if status == StatusCode::INTERNAL_SERVER_ERROR {
      error!("{:?}", self.0);
    }
    let body = json!({ "code": code, "message": format!("{:#}", self.0) });
    (status, Json(body)).into_response()
  }
}

impl<E: Into<Report>> From<E> for AppError {
  fn from(err: E) -> Self {
    Self(err.into())
  }
}
