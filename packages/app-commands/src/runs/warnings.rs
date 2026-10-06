use parking_lot::Mutex;
use treetime::progress::{LogLevel, LogSink, RunWarning};

pub struct WarningCollector<'a> {
  inner: &'a dyn LogSink,
  warnings: Mutex<Vec<RunWarning>>,
}

impl<'a> WarningCollector<'a> {
  pub fn new(inner: &'a dyn LogSink) -> Self {
    Self {
      inner,
      warnings: Mutex::new(vec![]),
    }
  }

  pub fn into_warnings(self) -> Vec<RunWarning> {
    self.warnings.into_inner()
  }
}

impl LogSink for WarningCollector<'_> {
  fn log(&self, level: LogLevel, message: &str) {
    self.inner.log(level, message);
  }

  fn log_enabled(&self, level: LogLevel) -> bool {
    self.inner.log_enabled(level)
  }

  fn warning(&self, warning: &RunWarning) {
    self.warnings.lock().push(warning.clone());
    self.inner.warning(warning);
  }
}
