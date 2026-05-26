use app_api::progress::{LogEvent, LogLevel, ProgressEvent, ProgressSink};
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use serde::Serialize;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

static CANCELLED: AtomicBool = AtomicBool::new(false);

pub fn cancel() {
  CANCELLED.store(true, Ordering::SeqCst);
}

pub fn reset_cancel() {
  CANCELLED.store(false, Ordering::SeqCst);
}

#[derive(Serialize)]
#[serde(tag = "type", content = "data")]
enum NapiEvent {
  #[serde(rename = "progress")]
  Progress(ProgressEvent),
  #[serde(rename = "log")]
  Log(LogEvent),
}

pub struct NapiProgressSink {
  tsfn: Arc<ThreadsafeFunction<String, ()>>,
}

impl NapiProgressSink {
  pub fn new(tsfn: Arc<ThreadsafeFunction<String, ()>>) -> Self {
    Self { tsfn }
  }

  fn send_event(&self, event: &NapiEvent) {
    if let Ok(json) = serde_json::to_string(&event) {
      self.tsfn.call(Ok(json), ThreadsafeFunctionCallMode::NonBlocking);
    }
  }
}

impl ProgressSink for NapiProgressSink {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    self.send_event(&NapiEvent::Progress(ProgressEvent {
      stage: stage.to_owned(),
      fraction,
      message: message.to_owned(),
    }));
  }

  fn log(&self, level: LogLevel, message: &str) {
    self.send_event(&NapiEvent::Log(LogEvent {
      level,
      message: message.to_owned(),
    }));
  }

  fn log_enabled(&self, _level: LogLevel) -> bool {
    true
  }

  fn is_cancelled(&self) -> bool {
    CANCELLED.load(Ordering::SeqCst)
  }
}
