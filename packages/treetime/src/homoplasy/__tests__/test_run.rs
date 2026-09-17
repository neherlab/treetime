#[cfg(test)]
mod tests {
  use crate::error::OperationError;
  use crate::homoplasy::pipeline::{self, HomoplasyInput, HomoplasyParams};

  #[test]
  fn test_run_returns_typed_not_implemented() {
    let err = pipeline::run(&HomoplasyParams, HomoplasyInput).unwrap_err();
    assert!(
      matches!(err, OperationError::NotImplemented("homoplasy")),
      "homoplasy run must return the typed not-implemented error, got: {err:?}"
    );
  }
}
