#[cfg(feature = "clap")]
use clap::ValueHint;
use smart_default::SmartDefault;
use std::path::PathBuf;

#[derive(Debug, Clone, SmartDefault)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ConfigArgs {
  /// Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line
  /// flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from
  /// the working directory for standard input; a relative path in a flag resolves from the working directory
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Config"))]
  config: Option<PathBuf>,
}
