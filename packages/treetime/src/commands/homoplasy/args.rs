use crate::commands::ancestral::args::TreetimeAncestralArgs;
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeHomoplasyArgs {
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub ancestral_args: TreetimeAncestralArgs,

  /// Number of constant sites not included in alignment
  #[cfg_attr(feature = "clap", clap(long = "const"))]
  pub constant_sites: Option<usize>,

  /// rescale branch lengths
  #[cfg_attr(feature = "clap", clap(long))]
  pub rescale: bool,

  /// generate a more detailed report
  #[cfg_attr(feature = "clap", clap(long))]
  pub detailed: Option<String>,

  /// TSV file containing DRM info. Columns headers: GENOMIC_POSITION, ALT_BASE, DRUG, GENE, SUBSTITUTION
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub drms: Option<PathBuf>,

  /// number of mutations/nodes that are printed to screen
  #[cfg_attr(feature = "clap", clap(long, short = 'n', default_value_t = 10))]
  #[default = 10]
  pub num_mut: usize,
}
