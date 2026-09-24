use app_commands::command::AppCommand;
use app_commands::job::{JobEvent, JobHandle, JobId, JobProgress, JobRegistry, JobStarted, TerminalEvent, run_job};
use eyre::{Report, WrapErr};
use serde_json::Value;
use std::sync::Arc;

pub fn start_job(
  registry: &Arc<JobRegistry>,
  job_id: &str,
  command: &str,
  config_json: &str,
) -> Result<PendingJob, Report> {
  let command: AppCommand = command
    .parse()
    .wrap_err_with(|| format!("When reading the command name `{command}`"))?;
  let config: Value = serde_json::from_str(config_json).wrap_err("When reading the command configuration")?;
  let handle = registry.register(JobId::parse(job_id)?)?;
  Ok(PendingJob {
    handle,
    command,
    config,
  })
}

pub struct PendingJob {
  handle: JobHandle,
  command: AppCommand,
  config: Value,
}

impl PendingJob {
  pub fn run<F: Fn(JobEvent) + Send + Sync>(self, emit: F) -> TerminalEvent {
    emit(JobEvent::Started(JobStarted {
      job_id: self.handle.job_id().clone(),
      command: self.command,
    }));
    let progress = JobProgress::new(emit);
    run_job(
      self.handle.job_id(),
      self.command,
      &self.config,
      &|_config: &mut Value| Ok(()),
      self.handle.token(),
      &progress,
    )
  }
}
