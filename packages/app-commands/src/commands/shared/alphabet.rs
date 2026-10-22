use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime::alphabet::alphabet::AlphabetName;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

/// Alphabet of the sequences, shared by every command that reads sequences.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlphabetArgs {
  /// Alphabet of the sequences
  ///
  /// `nuc` (default): nucleotides A, C, G and T, with the IUPAC ambiguity codes and `N` for an unknown base. `aa`: the
  /// 20 amino acids and the stop codon `*` as a 21st state, with `X` for an unknown residue. `aa-no-stop`: the 20
  /// amino acids without a state for the stop codon.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, hide_short_help = true, help_heading = "Input data")
  )]
  pub(crate) alphabet: Option<AlphabetNameCli>,
}

impl AlphabetArgs {
  pub fn alphabet_name(&self) -> Option<AlphabetName> {
    self.alphabet.map(Into::into)
  }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "AlphabetName")]
pub enum AlphabetNameCli {
  #[default]
  Nuc,
  Aa,
  #[cfg_attr(feature = "clap", value(name = "aa-no-stop"))]
  AaNoStop,
}

impl_from_mirror!(AlphabetNameCli => AlphabetName { Nuc, Aa, AaNoStop });
