use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use treetime::ancestral::params::MethodAncestral;

impl From<MethodAncestralCli> for MethodAncestral {
  fn from(method: MethodAncestralCli) -> Self {
    match method {
      MethodAncestralCli::Marginal => MethodAncestral::Marginal,
      MethodAncestralCli::Parsimony => MethodAncestral::Parsimony,
      MethodAncestralCli::Joint => MethodAncestral::Joint,
    }
  }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "MethodAncestral")]
pub enum MethodAncestralCli {
  #[default]
  Marginal,
  Parsimony,
  Joint,
}
