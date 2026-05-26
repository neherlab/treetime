use crate::commands::ancestral::args::TreetimeAncestralArgs;
use clap::{Parser, ValueHint};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Parser, Debug, SmartDefault, Serialize, Deserialize)]
pub struct TreetimeHomoplasyArgs {
  #[clap(flatten)]
  pub ancestral_args: TreetimeAncestralArgs,

  /// Number of constant sites not included in alignment
  #[clap(long = "const")]
  pub constant_sites: Option<usize>,

  /// rescale branch lengths
  #[clap(long)]
  pub rescale: bool,

  /// generate a more detailed report
  #[clap(long)]
  pub detailed: Option<String>,

  /// TSV file containing DRM info. Columns headers: GENOMIC_POSITION, ALT_BASE, DRUG, GENE, SUBSTITUTION
  #[clap(long)]
  #[clap(value_hint = ValueHint::FilePath)]
  pub drms: Option<PathBuf>,

  /// number of mutations/nodes that are printed to screen
  #[clap(long, short = 'n', default_value_t = 10)]
  #[default = 10]
  pub num_mut: usize,
}
