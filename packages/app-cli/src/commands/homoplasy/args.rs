use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeHomoplasyArgsRaw {
  #[cfg_attr(feature = "clap", clap(flatten))]
  pub ancestral_args: TreetimeAncestralArgsRaw,

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

/// Homoplasy arguments with required inputs proven present.
///
/// Produced from [`TreetimeHomoplasyArgsRaw`] by [`TryFrom`], which converts the embedded ancestral
/// args and so enforces ancestral's required `tree`.
#[derive(Debug, Clone)]
pub struct TreetimeHomoplasyArgs {
  pub ancestral_args: TreetimeAncestralArgs,
  pub constant_sites: Option<usize>,
  pub rescale: bool,
  pub detailed: Option<String>,
  pub drms: Option<PathBuf>,
  pub num_mut: usize,
}

impl TryFrom<TreetimeHomoplasyArgsRaw> for TreetimeHomoplasyArgs {
  type Error = Report;

  fn try_from(raw: TreetimeHomoplasyArgsRaw) -> Result<Self, Report> {
    Ok(Self {
      ancestral_args: TreetimeAncestralArgs::try_from(raw.ancestral_args)?,
      constant_sites: raw.constant_sites,
      rescale: raw.rescale,
      detailed: raw.detailed,
      drms: raw.drms,
      num_mut: raw.num_mut,
    })
  }
}
