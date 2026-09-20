use crate::cancel::CancelledError;
use derive_more::Display;
use eyre::Report;

#[derive(Debug, Display)]
pub enum OperationError {
  #[display("{_0}")]
  InvalidParams(Report),

  #[display("{_0}")]
  InvalidInput(Report),

  #[display("Operation cancelled")]
  Cancelled,

  #[display("{_0}")]
  InferenceFailed(Report),

  #[display("{_0}")]
  SinkFailed(Report),

  #[display("The {_0} operation is not yet implemented in v1")]
  NotImplemented(&'static str),
}

impl OperationError {
  pub fn into_report(self) -> Report {
    match self {
      Self::Cancelled => Report::new(CancelledError),
      Self::NotImplemented(operation) => Report::msg(format!("The {operation} operation is not yet implemented in v1")),
      Self::InvalidParams(report)
      | Self::InvalidInput(report)
      | Self::InferenceFailed(report)
      | Self::SinkFailed(report) => report,
    }
  }
}

impl std::error::Error for OperationError {}

impl From<Report> for OperationError {
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
