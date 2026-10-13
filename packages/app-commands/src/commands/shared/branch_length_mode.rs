use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime::optimize::params::BranchLengthMode;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "BranchLengthMode")]
pub enum BranchLengthModeCli {
  Input,
  #[default]
  Marginal,
}

impl_from_mirror!(BranchLengthModeCli => BranchLengthMode { Input, Marginal });
