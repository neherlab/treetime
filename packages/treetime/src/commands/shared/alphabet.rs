use crate::alphabet::alphabet::AlphabetName;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;

/// Alphabet selection shared by every command that reads sequences.
///
/// A single `--alphabet` flag replaces the earlier redundant pair of `--alphabet` and `--aa` (the
/// latter being a second way to request the amino-acid alphabet, which could disagree with
/// `--alphabet`). When the alphabet is not given, callers auto-detect it from sequence content and
/// fall back to the nucleotide alphabet when detection is ambiguous (see `detect_alphabet`).
///
/// The flag has no short form: `-a` is reserved for `--alignment`.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AlphabetArgs {
  /// Sequence alphabet
  ///
  /// When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when
  /// detection is ambiguous.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  pub alphabet: Option<AlphabetName>,
}
