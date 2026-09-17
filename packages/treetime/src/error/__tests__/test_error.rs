#[cfg(test)]
mod tests {
  use crate::cancel::CancelledError;
  use crate::error::OperationError;
  use eyre::{Report, WrapErr, eyre};
  use pretty_assertions::assert_eq;

  /// The messages of a report and its transitive causes, outermost first, as color_eyre renders them.
  fn cause_chain(report: &Report) -> Vec<String> {
    report.chain().map(ToString::to_string).collect()
  }

  #[test]
  fn test_error_into_report_preserves_cause_chain() {
    // A three-level report, the shape a pipeline `wrap_err` context produces.
    let original: Report = Err::<(), Report>(eyre!("inner cause"))
      .wrap_err("middle context")
      .wrap_err("outer context")
      .unwrap_err();

    let recovered = OperationError::InferenceFailed(original).into_report();

    // The full chain survives classification: no flattening to the top-level message.
    assert_eq!(
      vec!["outer context", "middle context", "inner cause"],
      cause_chain(&recovered)
    );
  }

  #[test]
  fn test_error_cancelled_into_report_rebuilds_cancelled_marker() {
    let recovered = OperationError::Cancelled.into_report();
    assert!(
      recovered.downcast_ref::<CancelledError>().is_some(),
      "recovered report must still downcast to CancelledError"
    );
    assert_eq!("Operation cancelled", recovered.to_string());
  }

  #[test]
  fn test_error_from_report_classifies_cancellation() {
    let classified = OperationError::from(Report::new(CancelledError));
    assert!(matches!(classified, OperationError::Cancelled));
  }

  #[test]
  fn test_error_from_report_defaults_to_inference_failed() {
    let classified = OperationError::from(eyre!("some numerical failure"));
    assert!(matches!(classified, OperationError::InferenceFailed(_)));
  }

  #[test]
  fn test_error_display_reports_carried_message() {
    let err = OperationError::InvalidParams(eyre!("bad flag {}", 42));
    assert_eq!("bad flag 42", err.to_string());
  }

  #[test]
  fn test_error_not_implemented_display_names_operation() {
    let err = OperationError::NotImplemented("homoplasy");
    assert_eq!("The homoplasy operation is not yet implemented in v1", err.to_string());
  }

  #[test]
  fn test_error_not_implemented_into_report_keeps_message() {
    let recovered = OperationError::NotImplemented("homoplasy").into_report();
    assert_eq!(
      "The homoplasy operation is not yet implemented in v1",
      recovered.to_string()
    );
  }
}
