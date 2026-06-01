#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::{Path, PathBuf};

/// Output destination shared by every command.
///
/// `--outdir` (short `-O`) names a directory into which the standard set of output files is written
/// under fixed names. Each command also exposes per-type `--output-*` flags that override the path of
/// a single output. An override path of `-` writes that output to standard output;
/// compression is chosen from the path extension by the underlying writer. This mirrors the Nextclade
/// scheme of explicit per-type output flags layered over a directory default: type is selected by the
/// flag, never inferred from the file extension.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct OutputArgs {
  /// Directory to write the standard set of output files to
  #[cfg_attr(feature = "clap", clap(long, short = 'O', value_hint = ValueHint::AnyPath))]
  pub outdir: PathBuf,
}

impl OutputArgs {
  /// Resolve a single output to its destination path.
  ///
  /// Returns the explicit per-type override when given, otherwise `outdir/<default_name>`.
  pub fn resolve(&self, override_path: Option<&Path>, default_name: &str) -> PathBuf {
    match override_path {
      Some(path) => path.to_path_buf(),
      None => self.outdir.join(default_name),
    }
  }
}
