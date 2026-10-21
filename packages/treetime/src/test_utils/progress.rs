use crate::progress::{LogEvent, LogLevel, LogSink, RunWarning};
use parking_lot::Mutex;

#[derive(Default)]
pub(crate) struct RecordingLog {
  events: Mutex<Vec<LogEvent>>,
  run_warnings: Mutex<Vec<RunWarning>>,
}

impl RecordingLog {
  pub(crate) fn warnings(&self) -> Vec<String> {
    self
      .events
      .lock()
      .iter()
      .filter(|event| event.level == LogLevel::Warn)
      .map(|event| event.message.clone())
      .collect()
  }

  pub(crate) fn run_warnings(&self) -> Vec<RunWarning> {
    self.run_warnings.lock().clone()
  }
}

impl LogSink for RecordingLog {
  fn log(&self, level: LogLevel, message: &str) {
    self.events.lock().push(LogEvent {
      level,
      message: message.to_owned(),
    });
  }

  fn log_enabled(&self, _level: LogLevel) -> bool {
    true
  }

  fn warning(&self, warning: &RunWarning) {
    self.log(LogLevel::Warn, &warning.message);
    self.run_warnings.lock().push(warning.clone());
  }
}
