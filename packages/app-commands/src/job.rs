use crate::command::{AppCommand, CommandOutcome};
use crate::json_float::JsonFloat;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::any::Any;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::atomic::{AtomicBool, Ordering};
use strum_macros::IntoStaticStr;
use treetime::cancel::{Cancel, CancelledError};
use treetime::progress::{LogEvent, LogLevel, LogSink, StageSink};
use treetime::timetree::convergence::metrics::IterationRecord;
use treetime_schema::ProgressEvent;
use treetime_utils::error::{ReportChain, panic_message};
use treetime_utils::make_error;

const JOB_ID_MAX_LEN: usize = 128;

/// Identifier of one command run, unique among the jobs of a process.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(try_from = "String")]
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

impl TryFrom<String> for JobId {
  type Error = Report;

  fn try_from(id: String) -> Result<Self, Report> {
    Self::parse(&id)
  }
}

/// Event of a running job, in the order the job emits them.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, IntoStaticStr)]
#[serde(tag = "type", content = "data", rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum JobEvent {
  /// The job was accepted; always the first event.
  Started(JobStarted),
  /// A stage of the computation began or advanced.
  Progress(ProgressEvent),
  /// A diagnostic message of the computation.
  Log(LogEvent),
  /// Convergence values of one timetree optimization iteration.
  Iteration(IterationEvent),
  /// The job ended; always the last event, exactly once per job.
  Terminal(TerminalEvent),
}

/// Convergence values of one timetree optimization iteration, as the tracelog records them, with the clock model the
/// iteration used.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct IterationEvent {
  /// Iteration number, starting at 1.
  pub iteration: usize,
  /// Number of ancestral sequence states that changed in this iteration.
  pub n_diff: usize,
  /// Number of nodes added by polytomy resolution in this iteration.
  pub n_resolved: usize,
  /// Largest change of a node time, in years.
  pub max_time_change: Option<JsonFloat>,
  /// Root-mean-square change of the node times, in years.
  pub rms_time_change: Option<JsonFloat>,
  /// Log likelihood of the sequences.
  pub log_lh_seq: Option<JsonFloat>,
  /// Log likelihood of the node positions under the clock model.
  pub log_lh_pos: Option<JsonFloat>,
  /// Log likelihood of the coalescent prior.
  pub log_lh_coal: Option<JsonFloat>,
  /// Sum of the available log likelihoods.
  pub log_lh_total: Option<JsonFloat>,
  /// Clock rate of the clock model the iteration used, in substitutions per site per year.
  pub clock_rate: JsonFloat,
  /// Squared correlation coefficient of the root-to-tip regression of that clock model; absent for a fixed rate.
  pub r_squared: Option<JsonFloat>,
}

impl From<&IterationRecord> for IterationEvent {
  fn from(record: &IterationRecord) -> Self {
    let metrics = &record.metrics;
    Self {
      iteration: record.iteration,
      n_diff: metrics.n_diff,
      n_resolved: metrics.n_resolved,
      max_time_change: metrics.max_time_change.map(JsonFloat),
      rms_time_change: metrics.rms_time_change.map(JsonFloat),
      log_lh_seq: metrics.log_lh_seq.map(JsonFloat::from),
      log_lh_pos: metrics.log_lh_pos.map(JsonFloat::from),
      log_lh_coal: metrics.log_lh_coal.map(JsonFloat::from),
      log_lh_total: metrics.log_lh_total.map(JsonFloat::from),
      clock_rate: JsonFloat(record.clock.clock_rate),
      r_squared: record.clock.r_squared.map(JsonFloat),
    }
  }
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
  /// The job stopped because the process that ran it stopped.
  Interrupted { job_id: JobId },
}

impl TerminalEvent {
  fn from_report(job_id: JobId, report: &Report, cancel: &dyn Cancel) -> Self {
    if report.downcast_ref::<CancelledError>().is_some() || cancel.is_cancelled() {
      return Self::Cancelled { job_id };
    }
    let ReportChain { message, causes } = ReportChain::of(report);
    Self::Error {
      job_id,
      message,
      causes,
    }
  }

  fn from_panic(job_id: JobId, payload: &(dyn Any + Send)) -> Self {
    let detail = panic_message(payload).unwrap_or_else(|| "no panic message".to_owned());
    Self::Error {
      job_id,
      message: format!("internal error: the computation panicked: {detail}"),
      causes: vec![],
    }
  }
}

pub fn run_job(
  job_id: &JobId,
  cancel: &dyn Cancel,
  work: impl FnOnce() -> Result<CommandOutcome, Report>,
) -> TerminalEvent {
  let outcome = catch_unwind(AssertUnwindSafe(work));
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

pub struct JobProgress<F: Fn(JobEvent) + Send + Sync> {
  emit: F,
}

impl<F: Fn(JobEvent) + Send + Sync> JobProgress<F> {
  pub const fn new(emit: F) -> Self {
    Self { emit }
  }
}

impl<F: Fn(JobEvent) + Send + Sync> StageSink for JobProgress<F> {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    (self.emit)(JobEvent::Progress(ProgressEvent {
      stage: stage.to_owned(),
      fraction,
      message: message.to_owned(),
    }));
  }

  fn iteration(&self, record: &IterationRecord) {
    (self.emit)(JobEvent::Iteration(IterationEvent::from(record)));
  }
}

impl<F: Fn(JobEvent) + Send + Sync> LogSink for JobProgress<F> {
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
