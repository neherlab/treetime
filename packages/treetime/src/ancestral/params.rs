use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use strum_macros::VariantNames;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, VariantNames, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[strum(serialize_all = "kebab-case")]
#[derive(Default)]
#[serde(rename_all = "kebab-case")]
pub enum MethodAncestral {
  #[default]
  Marginal,
  Parsimony,
  Joint,
}
