use crate::commands::shared::documented_default::documented_default;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub const REPRESENTATION_DEFAULT: RepresentationCli = RepresentationCli::Auto;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
#[cfg_attr(
  feature = "clap",
  command(mut_arg("representation", |arg| documented_default(&REPRESENTATION_DEFAULT).help(arg)))
)]
pub struct RepresentationArgs {
  /// Representation of the sequences on the tree
  ///
  /// The dense representation stores the probabilities of every state at every alignment position of each node; the
  /// sparse representation stores only the variable positions. Dense is more accurate when branches are long and many
  /// sites change, but uses more memory.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_enum, help_heading = "Ancestral reconstruction")
  )]
  #[schemars(transform = documented_default(&REPRESENTATION_DEFAULT))]
  pub representation: Option<RepresentationCli>,
}

pub fn representation_dense(representation: Option<RepresentationCli>) -> Option<bool> {
  representation.unwrap_or(REPRESENTATION_DEFAULT).dense()
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "Representation")]
pub enum RepresentationCli {
  /// Let TreeTime choose the representation; it chooses sparse
  #[default]
  Auto,
  /// Probabilities of every state at every alignment position
  Dense,
  /// Only the variable positions
  Sparse,
}

impl RepresentationCli {
  pub fn dense(self) -> Option<bool> {
    match self {
      Self::Auto => None,
      Self::Dense => Some(true),
      Self::Sparse => Some(false),
    }
  }
}
