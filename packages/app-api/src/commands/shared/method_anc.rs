use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use treetime::ancestral::params::MethodAncestral;

// CLI mirror of core `MethodAncestral` carrying the clap `ValueEnum` derive. Core keeps
// `MethodAncestral` as a plain domain enum with no command-line knowledge; the ancestral, timetree,
// and clock commands share this adapter copy for their `--method-anc` value parsing and convert back
// with `From`. Variants, serde spellings, and the `schemars(rename)` schema name are kept identical
// to the core enum so `--help`, config parsing, and the generated schema stay byte-identical. The
// comment is non-doc on purpose: the core enum carries no doc, so a doc comment here would add a
// `description` to the schema and break schema parity.
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

impl From<MethodAncestralCli> for MethodAncestral {
  fn from(method: MethodAncestralCli) -> Self {
    match method {
      MethodAncestralCli::Marginal => MethodAncestral::Marginal,
      MethodAncestralCli::Parsimony => MethodAncestral::Parsimony,
      MethodAncestralCli::Joint => MethodAncestral::Joint,
    }
  }
}
