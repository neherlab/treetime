use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime::optimize::params::BranchLengthMode;

#[derive(
  Copy,
  Debug,
  Clone,
  PartialEq,
  Eq,
  PartialOrd,
  Ord,
  SmartDefault,
  Serialize,
  Deserialize,
  JsonSchema,
  deser::Serialize,
  deser::Deserialize,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "BranchLengthMode")]
pub enum BranchLengthModeCli {
  Input,
  #[default]
  Marginal,
}

impl From<BranchLengthModeCli> for BranchLengthMode {
  fn from(mode: BranchLengthModeCli) -> Self {
    match mode {
      BranchLengthModeCli::Input => BranchLengthMode::Input,
      BranchLengthModeCli::Marginal => BranchLengthMode::Marginal,
    }
  }
}
