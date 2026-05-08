use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, ValueEnum, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
pub enum BranchLengthMode {
  Input,
  #[default]
  Marginal,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, ValueEnum, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
pub enum RerootMode {
  #[default]
  LeastSquares,
  MinDev,
  Oldest,
  ClockFilter,
  Mrca,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Serialize, Deserialize)]
#[value(rename_all = "kebab-case")]
#[derive(Default)]
pub enum MethodAncestral {
  Joint,
  #[default]
  Marginal,
  Parsimony,
}
