use tokio::sync::mpsc;
use treetime::progress::{LogLevel, ProgressEvent, ProgressSink};

pub struct ChannelProgress {
  tx: mpsc::UnboundedSender<ProgressEvent>,
}

impl ChannelProgress {
  pub fn new(tx: mpsc::UnboundedSender<ProgressEvent>) -> Self {
    Self { tx }
  }
}

impl ProgressSink for ChannelProgress {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    drop(self.tx.send(ProgressEvent {
      stage: stage.to_owned(),
      fraction,
      message: message.to_owned(),
    }));
  }

  fn log(&self, _level: LogLevel, _message: &str) {}
}
