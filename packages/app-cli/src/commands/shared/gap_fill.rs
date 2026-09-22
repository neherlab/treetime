use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime::seq::gap_fill::GapFill;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "GapFill")]
pub enum GapFillCli {
  #[default]
  OnlyTerminal,
  All,
  None,
}

impl From<GapFillCli> for GapFill {
  fn from(mode: GapFillCli) -> Self {
    match mode {
      GapFillCli::OnlyTerminal => GapFill::OnlyTerminal,
      GapFillCli::All => GapFill::All,
      GapFillCli::None => GapFill::None,
    }
  }
}

/// Gap-handling policy shared by every command that reads sequences.
///
/// Extracted from the per-command duplication of `gap_fill` plus `keep_overhangs` plus
/// `effective_gap_fill()`. The deprecated `--keep-overhangs` flag is retained, hidden, and overrides
/// `--gap-fill` to `none` for backward compatibility with v0 invocations.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GapFillArgs {
  /// How to handle gap characters in input sequences
  ///
  /// 'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0).
  /// 'all': replace all gap characters with the ambiguous character.
  /// 'none': leave all gap characters unchanged.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, default_value_t = GapFillCli::default(), conflicts_with = "keep_overhangs")
  )]
  gap_fill: GapFillCli,

  /// Do not fill terminal gaps (deprecated: use --gap-fill=none)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  keep_overhangs: bool,
}

impl GapFillArgs {
  pub(crate) fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill.into()
    }
  }
}
