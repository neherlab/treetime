#[cfg(test)]
mod tests {
  use crate::guard::{guarded, guarded_json, to_napi};
  use app_commands::bridge::error::{ErrorCode, ErrorResponse};
  use app_commands::runs::errors::not_found;
  use eyre::WrapErr;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_guard_passes_a_result_through() {
    assert_eq!(Ok(o!("{\"a\":1}")), guarded_json(|| Ok(json!({ "a": 1 }))));
  }

  #[test]
  fn test_guard_turns_an_error_into_a_typed_error_response() {
    let result: Result<(), _> = guarded(|| Err(not_found("no run `a`")).wrap_err("When reading run `a`"));
    let expected = ErrorResponse {
      code: ErrorCode::NotFound,
      message: o!("When reading run `a`"),
      causes: vec_of_owned!["no run `a`"],
    };
    assert_eq!(Err(expected), result);
  }

  #[test]
  fn test_guard_turns_a_panic_into_an_internal_error() {
    let result: Result<(), _> = guarded(|| panic!("index 3 out of range"));
    let expected = ErrorResponse {
      code: ErrorCode::InternalError,
      message: o!("the back end stopped the operation after an internal error"),
      causes: vec_of_owned!["index 3 out of range"],
    };
    assert_eq!(Err(expected), result);
  }

  #[test]
  fn test_guard_carries_the_error_response_as_json_in_the_napi_error() {
    let response = ErrorResponse {
      code: ErrorCode::Conflict,
      message: o!("run `a` has already started"),
      causes: vec![],
    };
    let error = to_napi(&response);
    assert_eq!(
      (
        napi::Status::GenericFailure,
        o!(r#"{"code":"conflict","message":"run `a` has already started","causes":[]}"#)
      ),
      (error.status, error.reason)
    );
  }
}
