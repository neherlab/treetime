use crate::timetree::convergence::metrics::IterationRecord;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use strum_macros::Display;

pub trait StageSink: Send + Sync {
  fn report(&self, stage: &str, fraction: f64, message: &str);
  fn iteration(&self, _record: &IterationRecord) {}
}

pub trait LogSink: Send + Sync {
  fn log(&self, level: LogLevel, message: &str);
  fn log_enabled(&self, level: LogLevel) -> bool;
  fn warning(&self, warning: &RunWarning) {
    if self.log_enabled(LogLevel::Warn) {
      self.log(LogLevel::Warn, &warning.message);
    }
  }
}

pub struct NoopProgress;

impl StageSink for NoopProgress {
  fn report(&self, _stage: &str, _fraction: f64, _message: &str) {}
}

impl LogSink for NoopProgress {
  fn log(&self, _level: LogLevel, _message: &str) {}
  fn log_enabled(&self, _level: LogLevel) -> bool {
    false
  }
}

#[derive(Debug, Clone, JsonSchema, Serialize, Deserialize)]
pub struct LogEvent {
  pub level: LogLevel,
  pub message: String,
}

/// A problem that the run found and continued past, which makes its results less reliable.
#[derive(Debug, Clone, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct RunWarning {
  /// Kind of the problem.
  pub kind: RunWarningKind,
  /// The problem, as a sentence; the run log shows the same text.
  pub message: String,
  /// Every name the problem concerns, sorted; an edge is named by its child node. The message may show only some of them.
  pub names: Vec<String>,
}

/// Kind of a run warning.
#[derive(Debug, Clone, Copy, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum RunWarningKind {
  /// More than one node of the input tree has the same name.
  DuplicateNodeNames,
  /// More than one sequence of an alignment has the same name.
  DuplicateSequenceNames,
  /// More than one row of the metadata table has the same name.
  DuplicateMetadataNames,
  /// More than one partition entry has the same name.
  DuplicatePartitionNames,
  /// The branch-length optimizer did not converge on some edges: it reached its iteration limit, or it stopped at the upper end of its search interval while the likelihood still increased; their lengths may not maximize the likelihood.
  BranchLengthsNotConverged,
  /// Some final dates lie outside the input dates: a tip outside its own input date range, or an internal node after the latest input date. The clock alone dates these nodes, as it does for the tips that the clock filter flagged.
  DatesOutsideInputDates,
  /// The coordinates file of a mugration run has no coordinates for some values of the trait; the map of the Auspice tree leaves these locations out.
  LocationsWithoutCoordinates,
}

#[derive(Debug, Clone, Copy, Display, PartialEq, Eq, PartialOrd, Ord, JsonSchema, Serialize, Deserialize)]
#[strum(serialize_all = "UPPERCASE")]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum LogLevel {
  Trace,
  Debug,
  Info,
  Warn,
  Error,
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
