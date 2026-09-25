#[cfg(test)]
mod tests {
  use crate::bridge::error::{ErrorCode, ErrorResponse};
  use crate::runs::errors::{UploadTooLarge, conflict, invalid, not_found};
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::panic::catch_unwind;
  use treetime_utils::{make_report, o, vec_of_owned};

  #[rustfmt::skip]
  #[rstest]
  #[case::not_found(not_found("no run `a`"),                                   ErrorCode::NotFound)]
  #[case::conflict( conflict("run `a` has already started"),                   ErrorCode::Conflict)]
  #[case::invalid(  invalid("bad title"),                                      ErrorCode::InvalidRequest)]
  #[case::upload(   Report::new(UploadTooLarge { limit: o!("1 B"), name: o!("x") }), ErrorCode::UploadTooLarge)]
  #[case::serde(    Report::new(serde_json::from_str::<u8>("x").unwrap_err()), ErrorCode::InvalidRequest)]
  #[case::internal( make_report!("disk full"),                                 ErrorCode::InternalError)]
  #[trace]
  fn test_error_code_classifies_the_error_types_of_the_run_layer(#[case] report: Report, #[case] expected: ErrorCode) {
    assert_eq!(expected, ErrorCode::of(&report));
  }

  #[test]
  fn test_error_code_finds_a_classified_error_under_added_context() {
    let report = not_found("no run `a`").wrap_err("When reading the run");
    assert_eq!(ErrorCode::NotFound, ErrorCode::of(&report));
  }

  #[test]
  fn test_error_response_splits_the_message_from_its_causes() {
    let report = Err::<(), _>(invalid("bad title"))
      .wrap_err("When renaming the run")
      .wrap_err("When updating run `a`")
      .unwrap_err();
    let expected = ErrorResponse {
      code: ErrorCode::InvalidRequest,
      message: o!("When updating run `a`"),
      causes: vec_of_owned!["When renaming the run", "bad title"],
    };
    assert_eq!(expected, ErrorResponse::from_report(&report));
  }

  #[test]
  fn test_error_response_serializes_the_code_in_snake_case() {
    let response = ErrorResponse {
      code: ErrorCode::UploadTooLarge,
      message: "too large".to_owned(),
      causes: vec![],
    };
    assert_eq!(
      r#"{"code":"upload_too_large","message":"too large","causes":[]}"#,
      serde_json::to_string(&response).unwrap()
    );
  }

  #[test]
  fn test_error_response_from_panic_keeps_the_panic_message() {
    let payload = catch_unwind(|| panic!("index {} out of range", 3)).unwrap_err();
    let expected = ErrorResponse {
      code: ErrorCode::InternalError,
      message: "the back end stopped the operation after an internal error".to_owned(),
      causes: vec!["index 3 out of range".to_owned()],
    };
    assert_eq!(expected, ErrorResponse::from_panic(payload.as_ref()));
  }

  #[test]
  fn test_error_response_from_panic_with_a_static_message() {
    let payload = catch_unwind(|| panic!("static message")).unwrap_err();
    assert_eq!(
      vec!["static message".to_owned()],
      ErrorResponse::from_panic(payload.as_ref()).causes
    );
  }
}
