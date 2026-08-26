#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Command configuration loaded from a file.
///
/// Every command flattens this so `--config` supplies the command's full configuration object,
/// serialized as JSON or YAML (optionally compressed: `.gz`, `.bz2`, `.xz`, `.zst`). The object has
/// the same shape the command serializes to, so a config dumped from one run can be replayed.
///
/// Precedence, highest to lowest: an explicit command-line flag, then the config file, then the
/// built-in default. The field itself is `#[serde(skip)]` so it never appears inside a config object
/// and a dumped configuration replays cleanly.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ConfigArgs {
  /// Load command configuration from a JSON or YAML file (optionally compressed).
  ///
  /// The file holds this command's configuration object, in the same shape the command serializes
  /// to. Format is chosen from the extension after any compression suffix (`.yaml`/`.yml` select
  /// YAML, otherwise JSON). Use `-` to read from stdin.
  ///
  /// Explicit command-line flags override values from the file; the file overrides defaults.
  #[serde(skip)]
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Config"))]
  pub config: Option<PathBuf>,
}
