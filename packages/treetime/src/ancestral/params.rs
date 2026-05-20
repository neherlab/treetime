use clap::ValueEnum;
use serde::{Deserialize, Serialize};

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
#[derive(Default)]
pub enum MethodAncestral {
  Joint,
  #[default]
  Marginal,
  Parsimony,
}
