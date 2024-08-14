use crate::io::json::{json_write_str, JsonPretty};
use crate::make_error;
use crate::utils::string::quote;
use clap::ArgEnum;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use indexmap::{IndexMap, IndexSet};
use itertools::{chain, Itertools};
use maplit::btreemap;
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use strum_macros::Display;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault, Display)]
#[clap(rename = "kebab-case")]
pub enum AlphabetName {
  #[default]
  Nuc,
  Aa,
}

pub type ProfileMap = IndexMap<char, Array1<f64>>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Alphabet {
  all: IndexSet<char>,
  canonical: BTreeSet<char>,
  ambiguous: BTreeMap<char, Vec<char>>,
  unknown: char,
  gap: char,
  profile_map: ProfileMap,
}

impl Alphabet {
  /// Create one of the pre-defined alphabets
  pub fn new(name: AlphabetName) -> Result<Self, Report> {
    match name {
      AlphabetName::Nuc => Self::with_config(&AlphabetConfig {
        canonical: vec!['A', 'C', 'G', 'T'],
        ambiguous: btreemap! {
          'R' => vec!['A', 'G'],
          'Y' => vec!['C', 'T'],
          'S' => vec!['C', 'G'],
          'W' => vec!['A', 'T'],
          'K' => vec!['G', 'T'],
          'M' => vec!['A', 'C'],
          'D' => vec!['A', 'G', 'T'],
          'H' => vec!['A', 'C', 'T'],
          'B' => vec!['C', 'G', 'T'],
          'V' => vec!['A', 'C', 'G'],
        },
        unknown: 'N',
        gap: '-',
      }),
      AlphabetName::Aa => Self::with_config(&AlphabetConfig {
        canonical: vec![
          'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
        ],
        ambiguous: btreemap! {
          'B' => vec!['N', 'D'],
          'Z' => vec!['Q', 'E'],
          'J' => vec!['L', 'I'],
        },
        unknown: 'X',
        gap: '-',
      }),
    }
  }

  /// Create custom alphabet from a given config
  pub fn with_config(cfg: &AlphabetConfig) -> Result<Self, Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      gap,
    } = cfg;

    let canonical: BTreeSet<char> = canonical.iter().copied().collect();
    if canonical.is_empty() {
      return make_error!("When creating alphabet: canonical set of characters is empty. This is not allowed.");
    }

    let ambiguous: BTreeMap<char, Vec<char>> = ambiguous.to_owned();

    let all: IndexSet<char> = chain!(
      canonical.iter().copied(),
      ambiguous.keys().copied(),
      [*unknown, *gap].into_iter(),
    )
    .collect();

    let profile_map = cfg.create_profile_map()?;

    Ok(Self {
      all,
      canonical,
      ambiguous,
      unknown: *unknown,
      gap: *gap,
      profile_map,
    })
  }

  #[inline]
  pub fn get_profile(&self, c: char) -> &Array1<f64> {
    self
      .profile_map
      .get(&c)
      .ok_or_else(|| {
        format!(
          "When accessing profile map: Unknown character: '{c}'. Known characters: {}",
          self.profile_map.keys().join(", ")
        )
      })
      .unwrap()
  }

  pub fn get_code(&self, profile: &Array1<f64>) -> char {
    // TODO(perf): this mapping needs to be precomputed
    self
      .profile_map
      .iter()
      .find_map(|(&c, p)| (p == profile).then_some(c))
      .ok_or_else(|| format!("When accessing profile map: Unknown profile: '{profile}'"))
      .unwrap()
  }

  pub fn sequence_to_indices<'a>(&'a self, chars: impl Iterator<Item = char> + 'a) -> impl Iterator<Item = usize> + 'a {
    chars.map(|c| self.index(c))
  }

  pub fn indices_to_sequence<'a>(
    &'a self,
    indices: impl Iterator<Item = usize> + 'a,
  ) -> impl Iterator<Item = char> + 'a {
    indices.map(|i| self.char(i))
  }

  /// All existing characters (including 'unknown' and 'gap')
  pub fn chars(&self) -> impl Iterator<Item = char> + '_ {
    self.all.iter().copied()
  }

  /// Get char by index (indexed in the same order as given by `.chars()`)
  pub fn char(&self, index: usize) -> char {
    self.all[index]
  }

  /// Get index of a character (indexed in the same order as given by `.chars()`)
  pub fn index(&self, c: char) -> usize {
    self.all.get_index_of(&c).unwrap()
  }

  /// Check if character is in alphabet (including 'unknown' and 'gap')
  pub fn is_char(&self, c: char) -> bool {
    self.all.contains(&c)
  }

  pub fn n_chars(&self) -> usize {
    self.all.len()
  }

  /// Canonical (unambiguous) characters (e.g. 'A', 'C', 'G', 'T' in nuc alphabet)
  pub fn canonical(&self) -> impl Iterator<Item = char> + '_ {
    self.canonical.iter().copied()
  }

  /// Check is character is canonical
  pub fn is_canonical(&self, c: char) -> bool {
    self.canonical().contains(&c)
  }

  pub fn n_canonical(&self) -> usize {
    self.canonical.len()
  }

  /// Ambiguous characters (e.g. 'R', 'S' etc. in nuc alphabet)
  pub fn ambiguous(&self) -> impl Iterator<Item = char> + '_ {
    self.ambiguous.keys().copied()
  }

  /// Check if character is ambiguous (e.g. 'R', 'S' etc. in nuc alphabet)
  pub fn is_ambiguous(&self, c: char) -> bool {
    self.ambiguous().contains(&c)
  }

  pub fn n_ambiguous(&self) -> usize {
    self.ambiguous.len()
  }

  /// Determined characters: canonical or ambiguous
  pub fn determined(&self) -> impl Iterator<Item = char> + '_ {
    chain!(self.canonical(), self.ambiguous())
  }

  pub fn is_determined(&self, c: char) -> bool {
    self.determined().contains(&c)
  }

  pub fn n_determined(&self) -> usize {
    self.determined().count()
  }

  /// Undetermined characters: gap or unknown
  pub fn undetermined(&self) -> impl Iterator<Item = char> + '_ {
    [self.gap(), self.unknown()].into_iter()
  }

  pub fn is_undetermined(&self, c: char) -> bool {
    self.undetermined().contains(&c)
  }

  pub fn n_undetermined(&self) -> usize {
    self.undetermined().count()
  }

  /// Get 'unknown' character
  pub fn unknown(&self) -> char {
    self.unknown
  }

  /// Check if character is an 'unknown' character
  pub fn is_unknown(&self, c: char) -> bool {
    c == self.unknown()
  }

  /// Get 'gap' character
  pub fn gap(&self) -> char {
    self.gap
  }

  /// Check if character is a gap
  pub fn is_gap(&self, c: char) -> bool {
    c == self.gap()
  }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct AlphabetConfig {
  pub canonical: Vec<char>,
  pub ambiguous: BTreeMap<char, Vec<char>>,
  pub unknown: char,
  pub gap: char,
}

impl AlphabetConfig {
  pub fn create_profile_map(&self) -> Result<ProfileMap, Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      ..
    } = self;

    self
      .validate()
      .wrap_err("When validating alphabet config")
      .with_section(|| {
        json_write_str(&self, JsonPretty(true))
          .unwrap()
          .header("Alphabet config")
      })?;

    let eye = Array2::<f64>::eye(canonical.len());

    // Add canonical to profile map
    let mut profile_map: ProfileMap = canonical
      .iter()
      .zip(eye.rows())
      .map(|(s, x)| (*s, x.to_owned()))
      .collect();

    // Add unknown to profile map
    profile_map.insert(*unknown, Array1::<f64>::ones(canonical.len()));

    // Add ambiguous to profile map
    ambiguous.iter().for_each(|(&key, values)| {
      let profile = canonical
        .iter()
        .enumerate()
        .map(|(i, c)| if values.contains(c) { 1.0 } else { 0.0 })
        .collect::<Array1<f64>>();
      profile_map.insert(key, profile);
    });

    Ok(profile_map)
  }

  pub fn validate(&self) -> Result<(), Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      gap,
    } = self;

    {
      let canonical_dupes = canonical.iter().duplicates().collect_vec();
      if !canonical_dupes.is_empty() {
        let msg = canonical_dupes.into_iter().join(", ");
        return make_error!("Canonical set contains duplicates: {msg}");
      }
    }

    {
      let ambiguous_dupes = ambiguous.keys().duplicates().collect_vec();
      if !ambiguous_dupes.is_empty() {
        let msg = ambiguous_dupes.into_iter().join(", ");
        return make_error!("Ambiguous set contains duplicates: {msg}");
      }
    }

    let canonical: BTreeSet<_> = canonical.iter().copied().collect();
    let ambiguous_keys: BTreeSet<_> = ambiguous.keys().copied().collect();
    let ambiguous_set_map: BTreeMap<char, BTreeSet<char>> = ambiguous
      .iter()
      .map(|(key, vals)| (*key, vals.iter().copied().collect()))
      .collect();
    {
      let canonical_inter_ambig: BTreeSet<_> = canonical.intersection(&ambiguous_keys).copied().collect();
      if !canonical_inter_ambig.is_empty() {
        let msg = canonical_inter_ambig.into_iter().join(", ");
        return make_error!("Canonical and ambiguous sets must be disjoint, but these characters are shared: {msg}");
      }
    }

    if canonical.contains(gap) {
      let msg = canonical.iter().map(quote).join(", ");
      return make_error!("Canonical set contains 'gap' character: {msg}");
    }

    if canonical.contains(unknown) {
      let msg = canonical.iter().map(quote).join(", ");
      return make_error!("Canonical set contains 'unknown' character: {msg}");
    }

    if ambiguous.keys().contains(&gap) {
      let msg = ambiguous.keys().map(quote).join(", ");
      return make_error!("Ambiguous set contains 'gap' character: {msg}");
    }

    if ambiguous.keys().contains(&gap) {
      let msg = ambiguous.keys().map(quote).join(", ");
      return make_error!("Ambiguous set contains 'unknown' character: {msg}");
    }

    {
      let ambig_gaps = ambiguous_set_map
        .iter()
        .map(|(key, vals)| (key, vals.difference(&canonical).collect::<BTreeSet<_>>()))
        .filter(|(key, extra)| !extra.is_empty())
        .collect_vec();
      if !ambig_gaps.is_empty() {
        let msg = ambig_gaps
          .iter()
          .map(|(key, vals)| format!("\"{key}\" => {}", vals.iter().map(quote).join(", ")))
          .join("; ");
        return make_error!("Ambiguous set can only contain mapping(s) to canonical characters, but found: {msg}");
      }
    }

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use indoc::indoc;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_alphabet_sequence_to_indices() -> Result<(), Report> {
    let actual = Alphabet::new(AlphabetName::Nuc)?
      .sequence_to_indices(array!['A', 'G', 'T', 'G', '-', 'G', 'N', 'G', 'C'].into_iter())
      .collect_vec();
    let expected = vec![0, 2, 3, 2, 15, 2, 14, 2, 1];
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_alphabet_indices_to_sequence() -> Result<(), Report> {
    let actual = Alphabet::new(AlphabetName::Nuc)?
      .indices_to_sequence(array![0, 2, 3, 2, 15, 2, 14, 2, 1].into_iter())
      .collect_vec();
    let expected = vec!['A', 'G', 'T', 'G', '-', 'G', 'N', 'G', 'C'];
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_alphabet_nuc() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let actual = json_write_str(&alphabet, JsonPretty(true))?;
    let expected = indoc! {r#"
    {
      "all": [
        "A",
        "C",
        "G",
        "T",
        "B",
        "D",
        "H",
        "K",
        "M",
        "R",
        "S",
        "V",
        "W",
        "Y",
        "N",
        "-"
      ],
      "canonical": [
        "A",
        "C",
        "G",
        "T"
      ],
      "ambiguous": {
        "B": [
          "C",
          "G",
          "T"
        ],
        "D": [
          "A",
          "G",
          "T"
        ],
        "H": [
          "A",
          "C",
          "T"
        ],
        "K": [
          "G",
          "T"
        ],
        "M": [
          "A",
          "C"
        ],
        "R": [
          "A",
          "G"
        ],
        "S": [
          "C",
          "G"
        ],
        "V": [
          "A",
          "C",
          "G"
        ],
        "W": [
          "A",
          "T"
        ],
        "Y": [
          "C",
          "T"
        ]
      },
      "unknown": "N",
      "gap": "-",
      "profile_map": {
        "A": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "B": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            1.0,
            1.0,
            1.0
          ]
        },
        "C": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            1.0,
            0.0,
            0.0
          ]
        },
        "D": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            0.0,
            1.0,
            1.0
          ]
        },
        "G": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            0.0,
            1.0,
            0.0
          ]
        },
        "H": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            1.0,
            0.0,
            1.0
          ]
        },
        "K": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            0.0,
            1.0,
            1.0
          ]
        },
        "M": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            1.0,
            0.0,
            0.0
          ]
        },
        "N": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            1.0,
            1.0,
            1.0
          ]
        },
        "R": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            0.0,
            1.0,
            0.0
          ]
        },
        "S": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            1.0,
            1.0,
            0.0
          ]
        },
        "T": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            1.0
          ]
        },
        "V": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            1.0,
            1.0,
            0.0
          ]
        },
        "W": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            1.0,
            0.0,
            0.0,
            1.0
          ]
        },
        "Y": {
          "v": 1,
          "dim": [
            4
          ],
          "data": [
            0.0,
            1.0,
            0.0,
            1.0
          ]
        }
      }
    }"#};
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_alphabet_aa() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Aa)?;
    let actual = json_write_str(&alphabet, JsonPretty(true))?;
    let expected = indoc! {r#"
    {
      "all": [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        "B",
        "J",
        "Z",
        "X",
        "-"
      ],
      "canonical": [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y"
      ],
      "ambiguous": {
        "B": [
          "N",
          "D"
        ],
        "J": [
          "L",
          "I"
        ],
        "Z": [
          "Q",
          "E"
        ]
      },
      "unknown": "X",
      "gap": "-",
      "profile_map": {
        "A": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "B": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "C": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "D": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "E": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "F": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "G": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "H": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "I": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "J": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "K": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "L": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "M": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "N": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "P": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "Q": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "R": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "S": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "T": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0
          ]
        },
        "V": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0
          ]
        },
        "W": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0
          ]
        },
        "X": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0
          ]
        },
        "Y": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0
          ]
        },
        "Z": {
          "v": 1,
          "dim": [
            20
          ],
          "data": [
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        }
      }
    }"#};
    assert_eq!(expected, actual);
    Ok(())
  }
}
