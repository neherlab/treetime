use derive_more::{Display, Error};
use eyre::Report;

#[derive(Debug, Display, Error)]
#[display("the inputs of a run are limited to {limit} in total; `{name}` does not fit")]
pub struct UploadTooLarge {
  #[error(not(source))]
  pub limit: String,
  #[error(not(source))]
  pub name: String,
}

pub fn not_found(message: impl Into<String>) -> Report {
  Report::new(RunNotFound {
    message: message.into(),
  })
}

#[derive(Debug, Display, Error)]
#[display("{message}")]
pub struct RunNotFound {
  #[error(not(source))]
  pub message: String,
}

pub fn conflict(message: impl Into<String>) -> Report {
  Report::new(RunConflict {
    message: message.into(),
  })
}

#[derive(Debug, Display, Error)]
#[display("{message}")]
pub struct RunConflict {
  #[error(not(source))]
  pub message: String,
}

pub fn invalid(message: impl Into<String>) -> Report {
  Report::new(InvalidRunRequest {
    message: message.into(),
  })
}

#[derive(Debug, Display, Error)]
#[display("{message}")]
pub struct InvalidRunRequest {
  #[error(not(source))]
  pub message: String,
}
