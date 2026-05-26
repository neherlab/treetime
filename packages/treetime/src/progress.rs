use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use strum_macros::Display;

#[derive(Debug, Clone, Copy, Display, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
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
  fn log_enabled(&self, level: LogLevel) -> bool;
}

pub struct NoopProgress;

impl ProgressSink for NoopProgress {
  fn report(&self, _stage: &str, _fraction: f64, _message: &str) {}
  fn log(&self, _level: LogLevel, _message: &str) {}
  fn log_enabled(&self, _level: LogLevel) -> bool {
    false
  }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ProgressEvent {
  pub stage: String,
  pub fraction: f64,
  pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct LogEvent {
  pub level: LogLevel,
  pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ErrorResponse {
  pub code: String,
  pub message: String,
}

#[macro_export]
macro_rules! progress_log {
  ($sink:expr, $level:expr, $($arg:tt)*) => {
    if $sink.log_enabled($level) {
      $sink.log($level, &format!($($arg)*));
    }
  };
}

#[macro_export]
macro_rules! progress_error {
  ($sink:expr, $($arg:tt)*) => {
    $crate::progress_log!($sink, $crate::progress::LogLevel::Error, $($arg)*)
  };
}

#[macro_export]
macro_rules! progress_warn {
  ($sink:expr, $($arg:tt)*) => {
    $crate::progress_log!($sink, $crate::progress::LogLevel::Warn, $($arg)*)
  };
}

#[macro_export]
macro_rules! progress_info {
  ($sink:expr, $($arg:tt)*) => {
    $crate::progress_log!($sink, $crate::progress::LogLevel::Info, $($arg)*)
  };
}

#[macro_export]
macro_rules! progress_debug {
  ($sink:expr, $($arg:tt)*) => {
    $crate::progress_log!($sink, $crate::progress::LogLevel::Debug, $($arg)*)
  };
}

#[macro_export]
macro_rules! progress_trace {
  ($sink:expr, $($arg:tt)*) => {
    $crate::progress_log!($sink, $crate::progress::LogLevel::Trace, $($arg)*)
  };
}
