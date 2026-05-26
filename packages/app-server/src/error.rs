use axum::Json;
use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use eyre::Report;

pub struct AppError(Report);

impl IntoResponse for AppError {
  fn into_response(self) -> Response {
    let body = serde_json::json!({
      "code": "internal_error",
      "message": format!("{:?}", self.0),
    });
    (StatusCode::INTERNAL_SERVER_ERROR, Json(body)).into_response()
  }
}

impl<E: Into<Report>> From<E> for AppError {
  fn from(err: E) -> Self {
    Self(err.into())
  }
}
