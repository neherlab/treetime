use app_commands::bridge::error::ErrorResponse;
use eyre::Report;
use serde::Serialize;
use std::panic::{AssertUnwindSafe, catch_unwind};

pub fn guarded<T>(operation: impl FnOnce() -> Result<T, Report>) -> Result<T, ErrorResponse> {
  match catch_unwind(AssertUnwindSafe(operation)) {
    Ok(Ok(value)) => Ok(value),
    Ok(Err(report)) => Err(ErrorResponse::from_report(&report)),
    Err(payload) => Err(ErrorResponse::from_panic(payload.as_ref())),
  }
}

pub fn guarded_json<T: Serialize>(operation: impl FnOnce() -> Result<T, Report>) -> Result<String, ErrorResponse> {
  guarded(|| Ok(serde_json::to_string(&operation()?)?))
}

pub fn to_napi(response: &ErrorResponse) -> napi::Error {
  let reason = serde_json::to_string(response).unwrap_or_else(|err| format!("{}: {err}", response.message));
  napi::Error::new(napi::Status::GenericFailure, reason)
}
