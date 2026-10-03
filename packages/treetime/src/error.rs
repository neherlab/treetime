use crate::cancel::CancelledError;
use derive_more::Display;
use eyre::Report;

#[derive(Debug, Display)]
#[display("{_0}")]
pub struct InputError(String);

impl std::error::Error for InputError {}

pub(crate) fn input_error(message: impl Into<String>) -> Report {
  Report::new(InputError(message.into()))
}

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

  pub fn classify(report: Report) -> Self {
    if report.downcast_ref::<CancelledError>().is_some() {
      Self::Cancelled
    } else if report.downcast_ref::<InputError>().is_some() {
      Self::InvalidInput(report)
    } else {
      Self::InferenceFailed(report)
    }
  }

  fn report(&self) -> Option<&Report> {
    match self {
      Self::InvalidParams(report)
      | Self::InvalidInput(report)
      | Self::InferenceFailed(report)
      | Self::SinkFailed(report) => Some(report),
      Self::Cancelled | Self::NotImplemented(_) => None,
    }
  }
}

impl std::error::Error for OperationError {
  fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
    let error: &(dyn std::error::Error + 'static) = self.report()?.as_ref();
    error.source()
  }
}
