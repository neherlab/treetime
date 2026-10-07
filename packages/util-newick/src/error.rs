use std::fmt::{self, Display};

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("{message}")]
pub struct NewickWriteError {
  message: String,
}

impl NewickWriteError {
  pub(crate) fn new(message: impl Into<String>) -> Self {
    Self {
      message: message.into(),
    }
  }
}

impl From<fmt::Error> for NewickWriteError {
  fn from(error: fmt::Error) -> Self {
    Self::new(error.to_string())
  }
}

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error(
  "{text:?} is not a Newick dialect: expected <structure>,<annotations> with structure one of {structures} and annotations one of {annotations}"
)]
pub struct ParseDialectError {
  pub text: String,
  pub structures: String,
  pub annotations: String,
}

pub(crate) trait WriteContext<T> {
  fn context(self, context: impl Display) -> Result<T, NewickWriteError>;

  fn with_context<C: Display>(self, context: impl FnOnce() -> C) -> Result<T, NewickWriteError>;
}

impl<T, E: Display> WriteContext<T> for Result<T, E> {
  fn context(self, context: impl Display) -> Result<T, NewickWriteError> {
    self.map_err(|error| NewickWriteError::new(format!("{context}: {error}")))
  }

  fn with_context<C: Display>(self, context: impl FnOnce() -> C) -> Result<T, NewickWriteError> {
    self.map_err(|error| NewickWriteError::new(format!("{}: {error}", context())))
  }
}

macro_rules! write_error {
  ($($arg:tt)*) => {
    $crate::error::NewickWriteError::new(format!($($arg)*))
  };
}

pub(crate) use write_error;
