use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime::optimize::params::BranchLengthMode;

// CLI mirror of core `BranchLengthMode` carrying the clap `ValueEnum` derive. Core keeps
// `BranchLengthMode` as a plain domain enum; the timetree and clock commands share this adapter copy
// for their `--branch-length-mode` value parsing and convert back with `From`. Variants, serde
// spellings, and the `schemars(rename)` schema name are kept identical to the core enum so `--help`,
// config parsing, and the generated schema stay byte-identical. The comment is non-doc on purpose:
// the core enum carries no doc, so a doc comment here would add a `description` to the schema and
// break schema parity.
#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
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
