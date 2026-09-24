use utoipa::ToSchema;

#[derive(ToSchema)]
pub struct DatasetInfo {
  pub name: String,
  pub files: Vec<String>,
}

#[derive(ToSchema)]
pub struct ErrorResponse {
  pub code: String,
  pub message: String,
}

#[derive(ToSchema, serde::Serialize)]
pub struct CancelJobResponse {
  pub cancelled: bool,
}
