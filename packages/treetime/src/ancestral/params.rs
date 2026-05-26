use serde::{Deserialize, Serialize};

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
#[derive(Default)]
pub enum MethodAncestral {
  Joint,
  #[default]
  Marginal,
  Parsimony,
}
