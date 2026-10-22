use crate::job::JobId;
use crate::runs::lifecycle::StartedRun;
use eyre::{Report, WrapErr};
use std::thread;

pub trait RunExecutor: Send + Sync {
  fn execute(&self, started: StartedRun) -> Result<(), Report>;

  fn cancel(&self, id: &JobId);
}

pub struct ThreadExecutor;

impl RunExecutor for ThreadExecutor {
  #[cfg_attr(
    dylint_lib = "custom",
    expect(
      spawn_handle_dropped,
      reason = "the run thread catches a panic of the run and records every outcome in the run store, which the app reads"
    )
  )]
  fn execute(&self, started: StartedRun) -> Result<(), Report> {
    let id = started.record().id.clone();
    thread::Builder::new()
      .name(format!("run-{}", id.as_str()))
      .spawn(move || started.run())
      .wrap_err_with(|| format!("When starting a thread for run `{}`", id.as_str()))?;
    Ok(())
  }

  fn cancel(&self, _id: &JobId) {}
}
