use app_commands::bridge::error::ErrorResponse;
use eyre::Report;
use std::panic::{AssertUnwindSafe, catch_unwind};
use treetime_utils::io::json::{JsonPretty, json_write_str};

pub(crate) fn guarded<T>(operation: impl FnOnce() -> Result<T, Report>) -> Result<T, ErrorResponse> {
  match catch_unwind(AssertUnwindSafe(operation)) {
    Ok(Ok(value)) => Ok(value),
    Ok(Err(report)) => Err(ErrorResponse::from_report(&report)),
    Err(payload) => Err(ErrorResponse::from_panic(payload.as_ref())),
  }
}

pub(crate) fn to_napi(response: &ErrorResponse) -> napi::Error {
  let reason = json_write_str(response, JsonPretty(false)).unwrap_or_else(|err| format!("{}: {err}", response.message));
  napi::Error::new(napi::Status::GenericFailure, reason)
}
