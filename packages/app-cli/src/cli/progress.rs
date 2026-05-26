use indicatif::{ProgressBar, ProgressStyle};
use parking_lot::Mutex;
use treetime::progress::{LogLevel, ProgressSink};

pub struct BarProgress {
  bar: ProgressBar,
  min_level: LogLevel,
}

impl BarProgress {
  pub fn new(min_level: LogLevel) -> Self {
    let bar = ProgressBar::new(1000);
    bar.set_style(
      ProgressStyle::with_template("{spinner:.green} [{bar:30.cyan/dim}] {percent}% {msg}")
        .expect("progress style template")
        .progress_chars("=> "),
    );
    Self { bar, min_level }
  }
}

impl Drop for BarProgress {
  fn drop(&mut self) {
    self.bar.finish_and_clear();
  }
}

impl ProgressSink for BarProgress {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    self.bar.set_position((fraction * 1000.0) as u64);
    if message.is_empty() {
      self.bar.set_message(stage.to_owned());
    } else {
      self.bar.set_message(format!("{stage}: {message}"));
    }
    if fraction >= 1.0 {
      self.bar.finish_and_clear();
    }
  }

  fn log(&self, level: LogLevel, message: &str) {
    if self.log_enabled(level) {
      self.bar.suspend(|| {
        eprintln!("[{level}] {message}");
      });
    }
  }

  fn log_enabled(&self, level: LogLevel) -> bool {
    level >= self.min_level
  }
}

pub struct TextProgress {
  min_level: LogLevel,
  last_stage: Mutex<String>,
}

impl TextProgress {
  pub fn new(min_level: LogLevel) -> Self {
    Self {
      min_level,
      last_stage: Mutex::new(String::new()),
    }
  }
}

impl ProgressSink for TextProgress {
  fn report(&self, stage: &str, _fraction: f64, _message: &str) {
    if self.log_enabled(LogLevel::Info) {
      let mut last = self.last_stage.lock();
      if *last != stage {
        *last = stage.to_owned();
        eprintln!("[INFO] {stage}");
      }
    }
  }

  fn log(&self, level: LogLevel, message: &str) {
    if self.log_enabled(level) {
      eprintln!("[{level}] {message}");
    }
  }

  fn log_enabled(&self, level: LogLevel) -> bool {
    level >= self.min_level
  }
}
