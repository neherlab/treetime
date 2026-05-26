//! Inspired by clap-verbosity-flag:
//! https://github.com/rust-cli/clap-verbosity-flag
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clap::{ArgAction, Args};
use log::LevelFilter;
use serde::Serialize;
use treetime::progress::LogLevel;

#[derive(Args, Debug, Clone)]
pub struct Verbosity {
  /// Set verbosity level of console output
  #[clap(long, global = true, value_parser = PossibleValuesParser::new(["off", "error", "warn", "info", "debug", "trace"])
      .map(|s| s.parse::<LevelFilter>().unwrap()))]
  #[clap(conflicts_with = "quiet", conflicts_with = "verbose", conflicts_with = "silent")]
  #[clap(default_value = "warn")]
  #[clap(display_order = 95)]
  pub verbosity: LevelFilter,

  /// Disable all console output. Same as `--verbosity=off`
  #[clap(long, global = true)]
  #[clap(conflicts_with = "quiet", conflicts_with = "verbose", conflicts_with = "verbosity")]
  #[clap(display_order = 96)]
  pub silent: bool,

  /// Make console output more verbose. Add multiple occurrences to increase verbosity further.
  #[clap(long, short = 'v', action = ArgAction::Count, global = true)]
  #[clap(conflicts_with = "quiet", conflicts_with = "verbosity", conflicts_with = "silent")]
  #[clap(display_order = 97)]
  pub verbose: u8,

  /// Make console output more quiet. Add multiple occurrences to make output even more quiet.
  #[clap(long, short = 'q', action = ArgAction::Count, global = true)]
  #[clap(conflicts_with = "verbose", conflicts_with = "verbosity")]
  #[clap(display_order = 98)]
  pub quiet: u8,

  /// Disable progress bar display
  #[clap(long, global = true)]
  #[clap(display_order = 99)]
  pub no_progress: bool,
}

impl Serialize for Verbosity {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    serializer.serialize_str(&self.get_filter_level().to_string().to_lowercase())
  }
}

impl Verbosity {
  pub const fn get_filter_level(&self) -> LevelFilter {
    if self.silent {
      LevelFilter::Off
    } else {
      let ilevel = level_to_int(self.verbosity);
      let ilevel = ilevel.saturating_add(self.verbose);
      let ilevel = ilevel.saturating_sub(self.quiet);
      level_from_int(ilevel)
    }
  }

  pub fn get_log_level(&self) -> Option<LogLevel> {
    match self.get_filter_level() {
      LevelFilter::Off => None,
      LevelFilter::Error => Some(LogLevel::Error),
      LevelFilter::Warn => Some(LogLevel::Warn),
      LevelFilter::Info => Some(LogLevel::Info),
      LevelFilter::Debug => Some(LogLevel::Debug),
      LevelFilter::Trace => Some(LogLevel::Trace),
    }
  }
}

const fn level_to_int(level: LevelFilter) -> u8 {
  match level {
    LevelFilter::Off => 0,
    LevelFilter::Error => 1,
    LevelFilter::Warn => 2,
    LevelFilter::Info => 3,
    LevelFilter::Debug => 4,
    LevelFilter::Trace => 5,
  }
}

const fn level_from_int(verbosity: u8) -> LevelFilter {
  match verbosity {
    0 => LevelFilter::Off,
    1 => LevelFilter::Error,
    2 => LevelFilter::Warn,
    3 => LevelFilter::Info,
    4 => LevelFilter::Debug,
    5.. => LevelFilter::Trace,
  }
}
