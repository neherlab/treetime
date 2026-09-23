use derive_more::{Display, Error};
use eyre::Report;

#[derive(Debug, Display, Error)]
#[display("Operation cancelled")]
pub struct CancelledError;

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
