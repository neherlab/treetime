use utoipa::ToSchema;

#[derive(ToSchema, serde::Serialize)]
pub struct ErrorResponse {
  pub code: String,
  pub message: String,
}

#[derive(ToSchema, serde::Serialize)]
pub struct CancelRunResponse {
  pub cancelled: bool,
}
