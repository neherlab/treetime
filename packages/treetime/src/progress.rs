use serde::{Deserialize, Serialize};
use strum_macros::Display;

#[derive(Debug, Clone, Copy, Display, PartialEq, Eq, Serialize, Deserialize)]
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
