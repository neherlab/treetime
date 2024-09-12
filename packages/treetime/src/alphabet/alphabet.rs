use crate::io::json::{json_write_str, JsonPretty};
use crate::make_error;
use crate::utils::string::quote;
use clap::ArgEnum;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use indexmap::{indexmap, IndexMap, IndexSet};
use itertools::{chain, Itertools};
use ndarray::{stack, Array1, Array2, Axis};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeSet;
use std::iter::once;
use strum_macros::Display;

pub const NON_CHAR: char = '.';
pub const VARIABLE_CHAR: char = '~';
pub const FILL_CHAR: char = ' ';

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
  canonical: IndexSet<char>,
  ambiguous: IndexMap<char, Vec<char>>,
  unknown: char,
  gap: char,
  treat_gap_as_unknown: bool,
  profile_map: ProfileMap,
}

impl Default for Alphabet {
  fn default() -> Self {
    Self::new(AlphabetName::Nuc, false).expect("Failed to create default alphabet")
  }
}

impl Alphabet {
  /// Create one of the pre-defined alphabets
  pub fn new(name: AlphabetName, treat_gap_as_unknown: bool) -> Result<Self, Report> {
    match name {
      AlphabetName::Nuc => Self::with_config(&AlphabetConfig {
        canonical: vec!['A', 'C', 'G', 'T'],
        ambiguous: indexmap! {
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
        treat_gap_as_unknown,
      }),
      AlphabetName::Aa => Self::with_config(&AlphabetConfig {
        canonical: vec![
          'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*',
        ],
        ambiguous: indexmap! {
          'B' => vec!['N', 'D'],
          'Z' => vec!['Q', 'E'],
          'J' => vec!['L', 'I'],
        },
        unknown: 'X',
        gap: '-',
        treat_gap_as_unknown,
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
      treat_gap_as_unknown,
    } = cfg;

    let canonical: IndexSet<char> = canonical.iter().copied().collect();
    if canonical.is_empty() {
      return make_error!("When creating alphabet: canonical set of characters is empty. This is not allowed.");
    }

    let ambiguous: IndexMap<char, Vec<char>> = ambiguous.to_owned();

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
      treat_gap_as_unknown: *treat_gap_as_unknown,
      profile_map,
    })
  }

  /// Resolve possible ambiguity of the given character to the set of canonical chars
  pub fn disambiguate(&self, c: char) -> BTreeSet<char> {
    // If unknown then could be any canonical (e.g. N => { A, C, G, T })
    if self.is_unknown(c) {
      self.canonical().collect()
    }
    // If ambiguous (e.g. R => { A, G })
    else if let Some(resolutions) = self.ambiguous.get(&c) {
      resolutions.iter().copied().collect()
    }
    // Otherwise it's not ambiguous and it's the char itself (incl. gap)
    else {
      once(c).collect()
    }
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

  #[allow(single_use_lifetimes)] // TODO: remove when anonymous lifetimes in `impl Trait` are stabilized
  pub fn seq2prof<'a>(&self, chars: impl IntoIterator<Item = &'a char>) -> Result<Array2<f64>, Report> {
    let prof = stack(
      Axis(0),
      &chars.into_iter().map(|&c| self.get_profile(c).view()).collect_vec(),
    )?;
    Ok(prof)
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
  pub fn contains(&self, c: char) -> bool {
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

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlphabetConfig {
  pub canonical: Vec<char>,
  pub ambiguous: IndexMap<char, Vec<char>>,
  pub unknown: char,
  pub gap: char,
  pub treat_gap_as_unknown: bool,
}

impl AlphabetConfig {
  pub fn create_profile_map(&self) -> Result<ProfileMap, Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      gap,
      treat_gap_as_unknown,
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

    if *treat_gap_as_unknown {
      // Add gap to profile map
      profile_map.insert(*gap, profile_map[unknown].clone());
    }

    Ok(profile_map)
  }

  pub fn validate(&self) -> Result<(), Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      gap,
      ..
    } = self;

    {
      let all = chain![
        self.canonical.iter().copied(),
        self.ambiguous.keys().copied(),
        self.ambiguous.values().flatten().copied(),
        once(self.unknown),
        once(self.gap)
      ]
      .collect_vec();

      for reserved in [NON_CHAR, VARIABLE_CHAR, FILL_CHAR] {
        if all.iter().any(|&c| c == reserved) {
          return make_error!("Alphabet contains reserved character: {reserved}");
        }
      }
    }

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

    let canonical: IndexSet<_> = canonical.iter().copied().collect();
    let ambiguous_keys: IndexSet<_> = ambiguous.keys().copied().collect();
    let ambiguous_set_map: IndexMap<char, IndexSet<char>> = ambiguous
      .iter()
      .map(|(key, vals)| (*key, vals.iter().copied().collect()))
      .collect();
    {
      let canonical_inter_ambig: IndexSet<_> = canonical.intersection(&ambiguous_keys).copied().collect();
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
        .map(|(key, vals)| (key, vals.difference(&canonical).collect::<IndexSet<_>>()))
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
  use maplit::btreeset;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_alphabet_sequence_to_indices() -> Result<(), Report> {
    let actual = Alphabet::new(AlphabetName::Nuc, false)?
      .sequence_to_indices(array!['A', 'G', 'T', 'G', '-', 'G', 'N', 'G', 'C'].into_iter())
      .collect_vec();
    let expected = vec![0, 2, 3, 2, 15, 2, 14, 2, 1];
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_alphabet_indices_to_sequence() -> Result<(), Report> {
    let actual = Alphabet::new(AlphabetName::Nuc, false)?
      .indices_to_sequence(array![0, 2, 3, 2, 15, 2, 14, 2, 1].into_iter())
      .collect_vec();
    let expected = vec!['A', 'G', 'T', 'G', '-', 'G', 'N', 'G', 'C'];
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_disambiguate() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    assert_eq!(btreeset! {'A', 'G'}, alphabet.disambiguate('R'));
    assert_eq!(btreeset! {'A', 'C', 'G', 'T'}, alphabet.disambiguate('N'));
    assert_eq!(btreeset! {'C'}, alphabet.disambiguate('C'));
    assert_eq!(btreeset! {alphabet.gap()}, alphabet.disambiguate(alphabet.gap()));
    Ok(())
  }

  #[test]
  fn test_alphabet_nuc() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let actual = json_write_str(&alphabet, JsonPretty(true))?;
    let expected = indoc! { /* language=json */ r#"{
      "all": [
        "A",
        "C",
        "G",
        "T",
        "R",
        "Y",
        "S",
        "W",
        "K",
        "M",
        "D",
        "H",
        "B",
        "V",
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
        "R": [
          "A",
          "G"
        ],
        "Y": [
          "C",
          "T"
        ],
        "S": [
          "C",
          "G"
        ],
        "W": [
          "A",
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
        "B": [
          "C",
          "G",
          "T"
        ],
        "V": [
          "A",
          "C",
          "G"
        ]
      },
      "unknown": "N",
      "gap": "-",
      "treat_gap_as_unknown": false,
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
        }
      }
    }"#};
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_alphabet_aa() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Aa, false)?;
    let actual = json_write_str(&alphabet, JsonPretty(true))?;
    let expected = indoc! { /* language=json */ r#"{
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
        "*",
        "B",
        "Z",
        "J",
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
        "Y",
        "*"
      ],
      "ambiguous": {
        "B": [
          "N",
          "D"
        ],
        "Z": [
          "Q",
          "E"
        ],
        "J": [
          "L",
          "I"
        ]
      },
      "unknown": "X",
      "gap": "-",
      "treat_gap_as_unknown": false,
      "profile_map": {
        "A": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "C": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "D": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "E": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "F": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "G": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "H": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "I": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "K": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "L": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "M": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "N": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "P": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "Q": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "R": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "S": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "T": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "V": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "W": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "Y": {
          "v": 1,
          "dim": [
            21
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
            1.0,
            0.0
          ]
        },
        "*": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            1.0
          ]
        },
        "X": {
          "v": 1,
          "dim": [
            21
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
            1.0,
            1.0
          ]
        },
        "B": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "Z": {
          "v": 1,
          "dim": [
            21
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
            0.0,
            0.0
          ]
        },
        "J": {
          "v": 1,
          "dim": [
            21
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
