use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime::ancestral::params::MethodAncestral;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "MethodAncestral")]
pub enum MethodAncestralCli {
  #[default]
  Marginal,
  Parsimony,
}

impl_from_mirror!(MethodAncestralCli => MethodAncestral { Marginal, Parsimony });
