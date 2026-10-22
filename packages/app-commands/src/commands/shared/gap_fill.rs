use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime::seq::gap_fill::GapFill;
use treetime_schema::schema_defaults;

/// Gap-handling policy shared by every command that reads sequences.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GapFillArgs {
  /// How to handle gap characters in input sequences
  ///
  /// 'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0).
  /// 'all': replace all gap characters with the ambiguous character.
  /// 'none': leave all gap characters unchanged.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, default_value_t = GapFillCli::default(), help_heading = "Ancestral reconstruction")
  )]
  gap_fill: GapFillCli,
}

impl GapFillArgs {
  pub fn gap_fill(&self) -> GapFill {
    self.gap_fill.into()
  }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "GapFill")]
pub enum GapFillCli {
  #[default]
  OnlyTerminal,
  All,
  None,
}

impl_from_mirror!(GapFillCli => GapFill { OnlyTerminal, All, None });
