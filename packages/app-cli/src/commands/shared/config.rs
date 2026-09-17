#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

/// Command configuration loaded from a file.
///
/// Every command flattens this so `--config` supplies the command's full configuration object. The
/// object has the same shape the command serializes to, so a config dumped from one run can be
/// replayed. The document is parsed as YAML; because YAML 1.2 is a superset of JSON, a JSON config
/// is accepted by the same parser without any extension-based format selection.
///
/// Precedence, highest to lowest: an explicit command-line flag, then the config file, then the
/// built-in default. The field itself is `#[serde(skip)]` so it never appears inside a config object
/// and a dumped configuration replays cleanly.
///
/// Boolean override is one-directional. A boolean flag has no command-line spelling for `false`, so a
/// boolean enabled in a config (for example `covariation: true`) cannot be turned off from the command
/// line; edit the config to disable it. An explicit flag can only raise a config `false` to `true`.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ConfigArgs {
  /// Load command configuration from a file (optionally compressed: `.gz`, `.bz2`, `.xz`, `.zst`).
  ///
  /// The file holds this command's configuration object, in the same shape the command serializes
  /// to. It is parsed as YAML, which also accepts JSON. Use `-` to read from stdin.
  ///
  /// Explicit command-line flags override values from the file; the file overrides defaults. A
  /// boolean enabled in the config cannot be disabled from the command line (a flag has no `false`
  /// spelling); edit the config instead.
  #[serde(skip)]
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Config"))]
  pub config: Option<PathBuf>,
}
