use crate::seq::gap_fill::GapFill;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;

/// Gap-handling policy shared by every command that reads sequences.
///
/// Extracted from the per-command duplication of `gap_fill` plus `keep_overhangs` plus
/// `effective_gap_fill()`. The deprecated `--keep-overhangs` flag is retained, hidden, and overrides
/// `--gap-fill` to `none` for backward compatibility with v0 invocations.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GapFillArgs {
  /// How to handle gap characters in input sequences
  ///
  /// 'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0).
  /// 'all': replace all gap characters with the ambiguous character.
  /// 'none': leave all gap characters unchanged.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, default_value_t = GapFill::default(), conflicts_with = "keep_overhangs")
  )]
  pub gap_fill: GapFill,

  /// Do not fill terminal gaps (deprecated: use --gap-fill=none)
  #[cfg_attr(feature = "clap", clap(long, hide = true))]
  pub keep_overhangs: bool,
}

impl GapFillArgs {
  /// Gap-fill mode after applying the deprecated `--keep-overhangs` override.
  pub fn effective_gap_fill(&self) -> GapFill {
    if self.keep_overhangs {
      GapFill::None
    } else {
      self.gap_fill
    }
  }
}
