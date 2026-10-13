use crate::error::{AppError, plain_error};
use crate::state::AppState;
use aide::axum::ApiRouter;
use app_backend::bridge::error::ErrorCode;
use axum::BoxError;
use axum::error_handling::HandleErrorLayer;
use axum::response::{IntoResponse, Response};
use eyre::eyre;
use std::sync::Arc;
use std::time::Duration;
use tower::ServiceBuilder;
use tower::timeout::error::Elapsed;

const REQUEST_TIMEOUT: Duration = Duration::from_secs(120);

pub(crate) fn with_request_timeout(routes: ApiRouter<Arc<AppState>>) -> ApiRouter<Arc<AppState>> {
  routes.layer(
    ServiceBuilder::new()
      .layer(HandleErrorLayer::new(timed_out))
      .timeout(REQUEST_TIMEOUT),
  )
}

async fn timed_out(error: BoxError) -> Response {
  if error.is::<Elapsed>() {
    plain_error(
      ErrorCode::Timeout,
      format!("the request took longer than {} seconds", REQUEST_TIMEOUT.as_secs()),
    )
  } else {
    AppError::from(eyre!(error)).into_response()
  }
}
