use crate::cancel::CancelledError;
use derive_more::Display;
use eyre::Report;

/// Typed failure class of a core operation.
///
/// Every core pipeline (`ancestral`, `optimize`, `prune`, `timetree`, `clock`, `mugration`) returns
/// this so a caller can react to the failure kind without matching on message text. Each fallible
/// variant carries the original [`Report`] unchanged, so the exact user-facing message and its cause
/// chain survive the classification and are recovered verbatim by [`OperationError::into_report`].
///
/// A caller that only needs the message and cause chain (the command-line path) must go through
/// [`OperationError::into_report`]: recovering the report directly preserves the eyre cause chain,
/// whereas letting eyre's blanket `From<E: Error>` wrap the value with `?` keeps only the top-level
/// message.
#[derive(Debug, Display)]
pub enum OperationError {
  /// A supplied parameter or flag combination is invalid.
  #[display("{_0}")]
  InvalidParams(Report),

  /// The input data (alignment, tree, dates, names, states) is missing or inconsistent.
  #[display("{_0}")]
  InvalidInput(Report),

  /// The operation observed a cancellation signal and stopped early.
  #[display("Operation cancelled")]
  Cancelled,

  /// A numerical or algorithmic stage of the inference failed.
  #[display("{_0}")]
  InferenceFailed(Report),

  /// A caller-supplied output sink failed while receiving a result.
  #[display("{_0}")]
  SinkFailed(Report),
}

impl OperationError {
  /// Recover the underlying [`Report`], preserving the exact message and cause chain.
  ///
  /// [`OperationError::Cancelled`] rebuilds a `CancelledError` report so downstream callers that
  /// downcast the returned report to [`CancelledError`] keep working and the `"Operation cancelled"`
  /// message is unchanged.
  pub fn into_report(self) -> Report {
    match self {
      Self::Cancelled => Report::new(CancelledError),
      Self::InvalidParams(report)
      | Self::InvalidInput(report)
      | Self::InferenceFailed(report)
      | Self::SinkFailed(report) => report,
    }
  }
}

// A plain error type with no `source`: the carried report holds the cause chain, recovered through
// `into_report` rather than the `Error::source` accessor.
impl std::error::Error for OperationError {}

impl From<Report> for OperationError {
  /// Classify a propagated report: a cancellation report becomes [`OperationError::Cancelled`],
  /// everything else defaults to [`OperationError::InferenceFailed`]. A pipeline constructs the other
  /// variants explicitly at the site where the failure class is known.
  fn from(report: Report) -> Self {
    if report.downcast_ref::<CancelledError>().is_some() {
      Self::Cancelled
    } else {
      Self::InferenceFailed(report)
    }
  }
}

#[cfg(test)]
mod __tests__;
