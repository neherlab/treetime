use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime::optimize::params::BranchLengthMode;

impl From<BranchLengthModeCli> for BranchLengthMode {
  fn from(mode: BranchLengthModeCli) -> Self {
    match mode {
      BranchLengthModeCli::Input => BranchLengthMode::Input,
      BranchLengthModeCli::Marginal => BranchLengthMode::Marginal,
    }
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "BranchLengthMode")]
pub(crate) enum BranchLengthModeCli {
  Input,
  #[default]
  Marginal,
}
