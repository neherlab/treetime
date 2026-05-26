pub use treetime::progress::{CancelledError, LogEvent, LogLevel, NoopProgress, ProgressEvent, ProgressSink};

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
