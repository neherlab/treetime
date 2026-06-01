use serde::{Deserialize, Serialize};
use strum_macros::VariantNames;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, VariantNames)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
#[strum(serialize_all = "kebab-case")]
#[derive(Default)]
pub enum MethodAncestral {
  #[default]
  Marginal,
  Parsimony,
  Joint,
}
