use eyre::{Report, WrapErr};
use indicatif::{ProgressBar, ProgressStyle};
use log::{Level, log};
use parking_lot::Mutex;
use treetime::progress::{LogLevel, LogSink, StageSink};

pub(crate) struct BarProgress {
  bar: ProgressBar,
  min_level: LogLevel,
}

impl BarProgress {
  pub(crate) fn new(min_level: LogLevel) -> Result<Self, Report> {
    let bar = ProgressBar::new(1000);
    bar.set_style(
      ProgressStyle::with_template("{spinner:.green} [{bar:30.cyan/dim}] {percent}% {msg}")
        .wrap_err("When parsing the progress bar template")?
        .progress_chars("=> "),
    );
    Ok(Self { bar, min_level })
  }
}

impl Drop for BarProgress {
  fn drop(&mut self) {
    self.bar.finish_and_clear();
  }
}

impl StageSink for BarProgress {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
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
}

impl LogSink for BarProgress {
  fn log(&self, level: LogLevel, message: &str) {
    if self.log_enabled(level) {
      self.bar.suspend(|| write_log(level, message));
    }
  }

  fn log_enabled(&self, level: LogLevel) -> bool {
    level >= self.min_level
  }
}

pub(crate) struct TextProgress {
  min_level: LogLevel,
  last_stage: Mutex<String>,
}

impl TextProgress {
  pub(crate) fn new(min_level: LogLevel) -> Self {
    Self {
      min_level,
      last_stage: Mutex::new(String::new()),
    }
  }
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    debug_remnants,
    reason = "the stage sink renders stage lines on stderr and has no error channel"
  )
)]
impl StageSink for TextProgress {
  fn report(&self, stage: &str, _fraction: f64, _message: &str) {
    if self.log_enabled(LogLevel::Info) {
      let mut last = self.last_stage.lock();
      if *last != stage {
        *last = stage.to_owned();
        eprintln!("[INFO] {stage}");
      }
    }
  }
}

impl LogSink for TextProgress {
  fn log(&self, level: LogLevel, message: &str) {
    if self.log_enabled(level) {
      write_log(level, message);
    }
  }

  fn log_enabled(&self, level: LogLevel) -> bool {
    level >= self.min_level
  }
}

fn write_log(level: LogLevel, message: &str) {
  let level = match level {
    LogLevel::Trace => Level::Trace,
    LogLevel::Debug => Level::Debug,
    LogLevel::Info => Level::Info,
    LogLevel::Warn => Level::Warn,
    LogLevel::Error => Level::Error,
  };
  log!(level, "{message}");
}
