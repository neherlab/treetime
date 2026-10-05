use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use smart_default::SmartDefault;
use treetime::alphabet::alphabet::AlphabetName;

/// Alphabet selection shared by every command that reads sequences.
///
/// When the alphabet is not given, callers auto-detect it from sequence content and fall back to the
/// nucleotide alphabet when detection is ambiguous (see `detect_alphabet`).
///
/// The flag has no short form: `-a` is reserved for `--alignment`.
#[skip_serializing_none]
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlphabetArgs {
  /// Sequence alphabet
  ///
  /// When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when
  /// detection is ambiguous.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Input data"))]
  alphabet: Option<AlphabetNameCli>,
}

impl AlphabetArgs {
  pub fn alphabet_name(&self) -> Option<AlphabetName> {
    self.alphabet.map(Into::into)
  }
}

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
