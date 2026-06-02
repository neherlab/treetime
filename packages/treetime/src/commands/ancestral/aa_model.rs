use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::get_gtr::GtrModelName;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use strum_macros::Display;

/// Amino-acid substitution model, mirroring the nucleotide `--model` but restricted to the values
/// that are sound over an amino-acid alphabet.
///
/// The default `infer` matches augur, which reconstructs amino acids with a JC69-seeded inferred
/// GTR over the stop-inclusive alphabet (`augur ancestral` calls `TreeAnc(..., gtr='JC69',
/// alphabet='aa')` with `infer_gtr=True`). Empirical matrices are opt-in.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Display, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum AaModelName {
  /// Infer an amino-acid GTR from the data over the stop-inclusive alphabet. Matches augur.
  #[default]
  Infer,
  /// Jones-Taylor-Thornton 1992 empirical 20-amino-acid model (no stop codon). Stop codons and any
  /// other out-of-alphabet characters in the input are mapped to the unknown state `X`.
  #[cfg_attr(feature = "clap", value(name = "jtt92"))]
  Jtt92,
}

impl AaModelName {
  /// Resolve to the concrete GTR model and the alphabet its reconstruction runs over.
  ///
  /// Inference is alphabet-agnostic, so `infer` runs over the 21-state stop-inclusive `Aa` alphabet
  /// and treats the stop codon as a real state. Empirical matrices have a fixed dimension, so
  /// `jtt92` runs over the 20-state `AaNoStop` alphabet and requires stops to be mapped to `X`.
  pub fn resolve(self) -> AaModelResolved {
    match self {
      AaModelName::Infer => AaModelResolved {
        gtr_model: GtrModelName::Infer,
        alphabet: AlphabetName::Aa,
      },
      AaModelName::Jtt92 => AaModelResolved {
        gtr_model: GtrModelName::Jtt92,
        alphabet: AlphabetName::AaNoStop,
      },
    }
  }
}

pub struct AaModelResolved {
  pub gtr_model: GtrModelName,
  pub alphabet: AlphabetName,
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_aa_model_default_is_infer() {
    assert_eq!(AaModelName::Infer, AaModelName::default());
  }

  #[test]
  fn test_aa_model_infer_resolves_to_stop_inclusive_inference() {
    let resolved = AaModelName::Infer.resolve();
    assert_eq!(GtrModelName::Infer, resolved.gtr_model);
    assert_eq!(AlphabetName::Aa, resolved.alphabet);
  }

  #[test]
  fn test_aa_model_jtt92_resolves_to_empirical_no_stop() {
    let resolved = AaModelName::Jtt92.resolve();
    assert_eq!(GtrModelName::Jtt92, resolved.gtr_model);
    assert_eq!(AlphabetName::AaNoStop, resolved.alphabet);
  }
}
