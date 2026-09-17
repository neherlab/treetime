use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime::alphabet::alphabet::AlphabetName;

// CLI mirror of core `AlphabetName` carrying the clap `ValueEnum` derive. Core keeps `AlphabetName` as
// a plain domain enum; this adapter copy owns the `--alphabet` value parsing and converts back with
// `From`. Variants, the `aa-no-stop` value spelling, serde spellings, and the `schemars(rename)`
// schema name are kept identical to the core enum so `--help`, config parsing, and the generated
// schema stay byte-identical. The comment is non-doc on purpose: the core enum carries no doc, so a
// doc comment here would add a `description` to the schema and break schema parity.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "AlphabetName")]
pub enum AlphabetNameCli {
  #[default]
  Nuc,
  Aa,
  #[cfg_attr(feature = "clap", value(name = "aa-no-stop"))]
  AaNoStop,
}

impl From<AlphabetNameCli> for AlphabetName {
  fn from(name: AlphabetNameCli) -> Self {
    match name {
      AlphabetNameCli::Nuc => AlphabetName::Nuc,
      AlphabetNameCli::Aa => AlphabetName::Aa,
      AlphabetNameCli::AaNoStop => AlphabetName::AaNoStop,
    }
  }
}

/// Alphabet selection shared by every command that reads sequences.
///
/// A single `--alphabet` flag replaces the earlier redundant pair of `--alphabet` and `--aa` (the
/// latter being a second way to request the amino-acid alphabet, which could disagree with
/// `--alphabet`). When the alphabet is not given, callers auto-detect it from sequence content and
/// fall back to the nucleotide alphabet when detection is ambiguous (see `detect_alphabet`).
///
/// The flag has no short form: `-a` is reserved for `--alignment`.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlphabetArgs {
  /// Sequence alphabet
  ///
  /// When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when
  /// detection is ambiguous.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  pub alphabet: Option<AlphabetNameCli>,
}

impl AlphabetArgs {
  /// Selected alphabet as the core domain value, or `None` for auto-detection.
  pub fn alphabet_name(&self) -> Option<AlphabetName> {
    self.alphabet.map(Into::into)
  }
}
