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
  /// How to handle the gap characters of the input sequences
  ///
  /// `only-terminal` (default, as in TreeTime 0.x): replace the leading and trailing gaps of each sequence with the
  /// unknown character of the alphabet (`N` or `X`). `all`: replace every gap with the unknown character. `none`: keep
  /// every gap.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_enum, default_value_t = GapFillCli::default(), help_heading = "Ancestral reconstruction")
  )]
  pub(crate) gap_fill: GapFillCli,
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
