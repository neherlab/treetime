<<<<<<< HEAD
<<<<<<< HEAD
pub use treetime::progress::{LogLevel, NoopProgress, ProgressSink};
=======
use strum_macros::Display;

#[derive(Debug, Clone, Copy, Display, PartialEq, Eq)]
#[strum(serialize_all = "UPPERCASE")]
pub enum LogLevel {
  Trace,
  Debug,
  Info,
  Warn,
  Error,
}

pub trait ProgressSink: Send + Sync {
  fn report(&self, stage: &str, fraction: f64, message: &str);
  fn log(&self, level: LogLevel, message: &str);
}

pub struct NoopProgress;

impl ProgressSink for NoopProgress {
  fn report(&self, _stage: &str, _fraction: f64, _message: &str) {}
  fn log(&self, _level: LogLevel, _message: &str) {}
}
>>>>>>> 7c23b04d (fix(api,napi,server): clean up Cargo deps, fix error contract, fix code style)
=======
pub use treetime::progress::{CancelledError, LogEvent, LogLevel, NoopProgress, ProgressEvent, ProgressSink};
>>>>>>> aad3402e (feat(app): Add progress reporting and cancellation support)

pub struct StderrProgress;

impl ProgressSink for StderrProgress {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    eprintln!("[{stage}] {:.0}% {message}", fraction * 100.0);
  }

  fn log(&self, level: LogLevel, message: &str) {
    eprintln!("[{level}] {message}");
  }

  fn log_enabled(&self, _level: LogLevel) -> bool {
    true
  }
}
