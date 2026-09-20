use eyre::Report;

#[derive(Debug)]
pub struct CancelledError;

impl std::fmt::Display for CancelledError {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    f.write_str("Operation cancelled")
  }
}

impl std::error::Error for CancelledError {}

pub trait Cancel: Send + Sync {
  fn is_cancelled(&self) -> bool;

  fn check(&self) -> Result<(), Report> {
    if self.is_cancelled() {
      Err(Report::new(CancelledError))
    } else {
      Ok(())
    }
  }
}

pub struct NoopCancel;

impl Cancel for NoopCancel {
  fn is_cancelled(&self) -> bool {
    false
  }
}
