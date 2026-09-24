use crate::command::{AppCommand, CommandOutcome};
use eyre::Report;
use parking_lot::Mutex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::any::Any;
use std::collections::BTreeMap;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use treetime::cancel::{Cancel, CancelledError};
use treetime::progress::{LogEvent, LogLevel, ProgressSink};
use treetime_schema::ProgressEvent;
use treetime_utils::make_error;

const JOB_ID_MAX_LEN: usize = 128;

/// Identifier of one command run, unique among the jobs of a process.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(transparent)]
pub struct JobId(String);

impl JobId {
  pub fn random() -> Self {
    Self(format!("{:032x}", rand::random::<u128>()))
  }

  pub fn parse(id: &str) -> Result<Self, Report> {
    let valid = !id.is_empty()
      && id.len() <= JOB_ID_MAX_LEN
      && id.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_');
    if valid {
      Ok(Self(id.to_owned()))
    } else {
      make_error!("invalid job id `{id}`: expected 1 to {JOB_ID_MAX_LEN} ASCII letters, digits, `-` or `_`")
    }
  }

  pub fn as_str(&self) -> &str {
    &self.0
  }
}

/// Event of a running job, in the order the job emits them.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "type", content = "data", rename_all = "kebab-case")]
pub enum JobEvent {
  /// The job was accepted; always the first event.
  Started(JobStarted),
  /// A stage of the computation began or advanced.
  Progress(ProgressEvent),
  /// A diagnostic message of the computation.
  Log(LogEvent),
  /// The job ended; always the last event, exactly once per job.
  Terminal(TerminalEvent),
}

/// Identity of an accepted job.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct JobStarted {
  pub job_id: JobId,
  pub command: AppCommand,
}

/// How a job ended.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "status", rename_all = "kebab-case")]
pub enum TerminalEvent {
  /// The command ran to completion.
  Ok { job_id: JobId, result: CommandOutcome },
  /// The command failed, or its configuration was rejected.
  Error {
    job_id: JobId,
    /// The error, as the CLI prints it.
    message: String,
    /// The errors that caused `message`, outermost first.
    causes: Vec<String>,
  },
  /// The job stopped because cancellation was requested.
  Cancelled { job_id: JobId },
}

impl TerminalEvent {
  fn from_report(job_id: JobId, report: &Report, cancel: &dyn Cancel) -> Self {
    if report.downcast_ref::<CancelledError>().is_some() || cancel.is_cancelled() {
      return Self::Cancelled { job_id };
    }
    Self::Error {
      job_id,
      message: report.to_string(),
      causes: report.chain().skip(1).map(ToString::to_string).collect(),
    }
  }

  fn from_panic(job_id: JobId, payload: &(dyn Any + Send)) -> Self {
    let detail = payload
      .downcast_ref::<&str>()
      .map(|message| (*message).to_owned())
      .or_else(|| payload.downcast_ref::<String>().cloned())
      .unwrap_or_else(|| "no panic message".to_owned());
    Self::Error {
      job_id,
      message: format!("internal error: the computation panicked: {detail}"),
      causes: vec![],
    }
  }
}

pub fn run_job(
  job_id: &JobId,
  command: AppCommand,
  config: &Value,
  prepare_config: &dyn Fn(&mut Value) -> Result<(), Report>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> TerminalEvent {
  let outcome = catch_unwind(AssertUnwindSafe(|| {
    let mut config = config.clone();
    prepare_config(&mut config)?;
    let prepared = command.prepare_value(&config)?;
    cancel.check()?;
    prepared.args.run(cancel, progress)
  }));
  match outcome {
    Ok(Ok(result)) => TerminalEvent::Ok {
      job_id: job_id.clone(),
      result,
    },
    Ok(Err(report)) => TerminalEvent::from_report(job_id.clone(), &report, cancel),
    Err(payload) => TerminalEvent::from_panic(job_id.clone(), payload.as_ref()),
  }
}

#[derive(Debug, Default)]
pub struct CancelToken(AtomicBool);

impl CancelToken {
  pub fn cancel(&self) {
    self.0.store(true, Ordering::SeqCst);
  }
}

impl Cancel for CancelToken {
  fn is_cancelled(&self) -> bool {
    self.0.load(Ordering::SeqCst)
  }
}

#[derive(Debug, Default)]
pub struct JobRegistry {
  jobs: Mutex<BTreeMap<JobId, Arc<CancelToken>>>,
}

impl JobRegistry {
  pub fn register(self: &Arc<Self>, job_id: JobId) -> Result<JobHandle, Report> {
    let token = Arc::new(CancelToken::default());
    let mut jobs = self.jobs.lock();
    if jobs.contains_key(&job_id) {
      return make_error!("a job with id `{}` is already running", job_id.as_str());
    }
    jobs.insert(job_id.clone(), Arc::clone(&token));
    Ok(JobHandle {
      registry: Arc::clone(self),
      job_id,
      token,
    })
  }

  pub fn cancel(&self, job_id: &JobId) -> bool {
    self.jobs.lock().get(job_id).is_some_and(|token| {
      token.cancel();
      true
    })
  }
}

#[derive(Debug)]
pub struct JobHandle {
  registry: Arc<JobRegistry>,
  job_id: JobId,
  token: Arc<CancelToken>,
}

impl JobHandle {
  pub const fn job_id(&self) -> &JobId {
    &self.job_id
  }

  pub fn token(&self) -> &CancelToken {
    &self.token
  }
}

impl Drop for JobHandle {
  fn drop(&mut self) {
    self.registry.jobs.lock().remove(&self.job_id);
  }
}

pub struct JobProgress<F: Fn(JobEvent) + Send + Sync> {
  emit: F,
}

impl<F: Fn(JobEvent) + Send + Sync> JobProgress<F> {
  pub const fn new(emit: F) -> Self {
    Self { emit }
  }
}

impl<F: Fn(JobEvent) + Send + Sync> ProgressSink for JobProgress<F> {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    (self.emit)(JobEvent::Progress(ProgressEvent {
      stage: stage.to_owned(),
      fraction,
      message: message.to_owned(),
    }));
  }

  fn log(&self, level: LogLevel, message: &str) {
    (self.emit)(JobEvent::Log(LogEvent {
      level,
      message: message.to_owned(),
    }));
  }

  fn log_enabled(&self, _level: LogLevel) -> bool {
    true
  }
}
