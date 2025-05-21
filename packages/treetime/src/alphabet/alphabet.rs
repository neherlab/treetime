use crate::io::json::{JsonPretty, json_write_str};
use crate::representation::bitset128::BitSet128;
use crate::representation::seq_char::AsciiChar;
use crate::representation::state_set::StateSet;
use crate::utils::string::quote;
use crate::{make_error, stateset, vec_u8};
use clap::ValueEnum;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use indexmap::{IndexMap, indexmap};
use itertools::{Itertools, chain};
use ndarray::{Array1, Array2, Axis, stack};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::borrow::Borrow;
use std::fmt::Display;
use std::iter::once;
use strum_macros::Display;

pub const NON_CHAR: AsciiChar = AsciiChar(b'.');
pub const VARIABLE_CHAR: AsciiChar = AsciiChar(b'~');
pub const FILL_CHAR: AsciiChar = AsciiChar(b' ');

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Display, ValueEnum)]
#[value(rename_all = "kebab-case")]
pub enum AlphabetName {
  #[default]
  Nuc,
  Aa,
}

pub type ProfileMap = IndexMap<AsciiChar, Array1<f64>>;
pub type StateSetMap = IndexMap<AsciiChar, StateSet>;
pub type CharToSet = IndexMap<AsciiChar, StateSet>;
pub type SetToChar = IndexMap<StateSet, AsciiChar>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Alphabet {
  all: StateSet,
  canonical: StateSet,
  ambiguous: IndexMap<AsciiChar, Vec<AsciiChar>>,
  ambiguous_keys: StateSet,
  determined: StateSet,
  undetermined: StateSet,
  unknown: AsciiChar,
  gap: AsciiChar,
  treat_gap_as_unknown: bool,
  profile_map: ProfileMap,

  #[serde(skip)]
  char_to_set: IndexMap<AsciiChar, StateSet>,
  #[serde(skip)]
  set_to_char: IndexMap<StateSet, AsciiChar>,

  #[serde(skip)]
  char_to_index: Vec<Option<usize>>,
  #[serde(skip)]
  index_to_char: Vec<AsciiChar>,
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
        canonical: vec_u8!['A', 'C', 'G', 'T'],
        ambiguous: indexmap! {
          b'R' => vec_u8!['A', 'G'],
          b'Y' => vec_u8!['C', 'T'],
          b'S' => vec_u8!['C', 'G'],
          b'W' => vec_u8!['A', 'T'],
          b'K' => vec_u8!['G', 'T'],
          b'M' => vec_u8!['A', 'C'],
          b'D' => vec_u8!['A', 'G', 'T'],
          b'H' => vec_u8!['A', 'C', 'T'],
          b'B' => vec_u8!['C', 'G', 'T'],
          b'V' => vec_u8!['A', 'C', 'G'],
        },
        unknown: b'N',
        gap: b'-',
        treat_gap_as_unknown,
      }),
      AlphabetName::Aa => Self::with_config(&AlphabetConfig {
        canonical: vec_u8![
          'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*',
        ],
        ambiguous: indexmap! {
          b'B' => vec_u8!['N', 'D'],
          b'Z' => vec_u8!['Q', 'E'],
          b'J' => vec_u8!['L', 'I'],
        },
        unknown: b'X',
        gap: b'-',
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

    let gap = AsciiChar::from(*gap);
    let unknown = AsciiChar::from(*unknown);

    let canonical = StateSet::from_iter(canonical);
    if canonical.is_empty() {
      return make_error!("When creating alphabet: canonical set of characters is empty. This is not allowed.");
    }

    let ambiguous: IndexMap<AsciiChar, Vec<AsciiChar>> = ambiguous
      .iter()
      .map(|(k, v)| (AsciiChar(*k), v.iter().copied().map(AsciiChar).collect()))
      .collect();
    let ambiguous_keys = ambiguous.keys().collect();

    let undetermined = stateset! {unknown, gap};
    let determined = StateSet::from_union([canonical, ambiguous_keys]);
    let all = StateSet::from_union([canonical, ambiguous_keys, undetermined]);

    let profile_map = cfg.create_profile_map()?;

    let mut char_to_index = vec![None; 128];
    let mut index_to_char = Vec::with_capacity(canonical.len());
    for (i, c) in canonical.iter().enumerate() {
      char_to_index[usize::from(c)] = Some(i);
      index_to_char.push(c);
    }

    let char_to_set = {
      let mut char_to_set: CharToSet = canonical.iter().map(|c| (c, StateSet::from_char(c))).collect();
      ambiguous.iter().for_each(|(key, chars)| {
        char_to_set.insert(*key, StateSet::from_iter(chars));
      });
      char_to_set.insert(gap, StateSet::from_char(gap));
      char_to_set.insert(unknown, StateSet::from_char(unknown));
      char_to_set
    };

    let set_to_char: SetToChar = char_to_set.iter().map(|(&c, &s)| (s, c)).collect();

    Ok(Self {
      all,
      char_to_index,
      index_to_char,
      canonical,
      ambiguous,
      ambiguous_keys,
      determined,
      undetermined,
      unknown,
      gap,
      treat_gap_as_unknown: *treat_gap_as_unknown,
      profile_map,
      char_to_set,
      set_to_char,
    })
  }

  #[inline]
  pub fn get_profile(&self, c: AsciiChar) -> &Array1<f64> {
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

  /// Create a profile vector given a set of characters
  pub fn construct_profile<I, T>(&self, chars: I) -> Result<Array1<f64>, Report>
  where
    I: IntoIterator<Item = T>,
    T: Borrow<AsciiChar> + Display,
  {
    let mut profile = Array1::<f64>::zeros(self.n_canonical());
    for c in chars {
      let chars = self.char_to_set(*c.borrow());
      for c in chars.iter() {
        let index = self.index(c);
        profile[index] = 1.0;
      }
    }
    Ok(profile)
  }

  pub fn get_code(&self, profile: &Array1<f64>) -> AsciiChar {
    // TODO(perf): this mapping needs to be precomputed
    self
      .profile_map
      .iter()
      .find_map(|(&c, p)| (p == profile).then_some(c))
      .ok_or_else(|| format!("When accessing profile map: Unknown profile: '{profile}'"))
      .unwrap()
  }

  #[allow(single_use_lifetimes)] // TODO: remove when anonymous lifetimes in `impl Trait` are stabilized
  pub fn seq2prof<'a>(&self, chars: impl IntoIterator<Item = &'a AsciiChar>) -> Result<Array2<f64>, Report> {
    let prof = stack(
      Axis(0),
      &chars.into_iter().map(|&c| self.get_profile(c).view()).collect_vec(),
    )?;
    Ok(prof)
  }

  pub fn set_to_char(&self, c: StateSet) -> AsciiChar {
    self.set_to_char[&c]
  }

  pub fn char_to_set(&self, c: impl Into<AsciiChar>) -> StateSet {
    self.char_to_set[&c.into()]
  }

  /// All existing characters (including 'unknown' and 'gap')
  pub fn chars(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    self.all.iter()
  }

  /// Get u8 by index (indexed in the same order as given by `.chars()`)
  pub fn char(&self, index: usize) -> AsciiChar {
    self.index_to_char[index]
  }

  /// Get index of a character (indexed in the same order as given by `.chars()`)
  pub fn index(&self, c: impl Into<usize>) -> usize {
    self.char_to_index[c.into()].unwrap()
  }

  /// Check if character is in alphabet (including 'unknown' and 'gap')
  pub fn contains(&self, c: AsciiChar) -> bool {
    self.all.contains(c)
  }

  pub fn n_chars(&self) -> usize {
    self.all.len()
  }

  /// Canonical (unambiguous) characters (e.g. 'A', 'C', 'G', 'T' in nuc alphabet)
  pub fn canonical(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    self.canonical.iter()
  }

  /// Check is character is canonical
  pub fn is_canonical(&self, c: AsciiChar) -> bool {
    self.canonical.contains(c)
  }

  pub fn n_canonical(&self) -> usize {
    self.canonical.len()
  }

  /// Ambiguous characters (e.g. 'R', 'S' etc. in nuc alphabet)
  pub fn ambiguous(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    self.ambiguous_keys.iter()
  }

  /// Check if character is ambiguous (e.g. 'R', 'S' etc. in nuc alphabet)
  pub fn is_ambiguous(&self, c: AsciiChar) -> bool {
    self.ambiguous_keys.contains(c)
  }

  pub fn n_ambiguous(&self) -> usize {
    self.ambiguous.len()
  }

  /// Determined characters: canonical or ambiguous
  pub fn determined(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    self.determined.iter()
  }

  pub fn is_determined(&self, c: AsciiChar) -> bool {
    self.determined.contains(c)
  }

  pub fn n_determined(&self) -> usize {
    self.determined.len()
  }

  /// Undetermined characters: gap or unknown
  pub fn undetermined(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    self.undetermined.iter()
  }

  pub fn is_undetermined(&self, c: AsciiChar) -> bool {
    self.undetermined.contains(c)
  }

  pub fn n_undetermined(&self) -> usize {
    self.undetermined.len()
  }

  /// Get 'unknown' character
  pub fn unknown(&self) -> AsciiChar {
    self.unknown
  }

  /// Check if character is an 'unknown' character
  pub fn is_unknown(&self, c: impl Into<AsciiChar>) -> bool {
    c.into() == self.unknown()
  }

  /// Get 'gap' character
  pub fn gap(&self) -> AsciiChar {
    self.gap
  }

  /// Check if character is a gap
  pub fn is_gap(&self, c: impl Into<AsciiChar>) -> bool {
    c.into() == self.gap()
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlphabetConfig {
  pub canonical: Vec<u8>,
  pub ambiguous: IndexMap<u8, Vec<u8>>,
  pub unknown: u8,
  pub gap: u8,
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

    let gap = AsciiChar::from(*gap);
    let unknown = AsciiChar::from(*unknown);

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
      .map(|(s, x)| (AsciiChar(*s), x.to_owned()))
      .collect();

    // Add unknown to profile map
    profile_map.insert(unknown, Array1::<f64>::ones(canonical.len()));

    // Add ambiguous to profile map
    ambiguous.iter().for_each(|(&key, values)| {
      let profile = canonical
        .iter()
        .enumerate()
        .map(|(i, c)| if values.contains(c) { 1.0 } else { 0.0 })
        .collect::<Array1<f64>>();
      profile_map.insert(AsciiChar(key), profile);
    });

    if *treat_gap_as_unknown {
      // Add gap to profile map
      profile_map.insert(gap, profile_map[&unknown].clone());
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
        if all.iter().any(|&c| c == u8::from(reserved)) {
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

    let canonical: StateSet = canonical.iter().copied().collect();
    let ambiguous_keys: StateSet = ambiguous.keys().copied().collect();
    let ambiguous_set_map: IndexMap<u8, StateSet> = ambiguous
      .iter()
      .map(|(key, vals)| (*key, vals.iter().copied().collect()))
      .collect();
    {
      let canonical_inter_ambig: StateSet = canonical.intersection(&ambiguous_keys);
      if !canonical_inter_ambig.is_empty() {
        let msg = canonical_inter_ambig.iter().join(", ");
        return make_error!("Canonical and ambiguous sets must be disjoint, but these characters are shared: {msg}");
      }
    }

    if canonical.contains(*gap) {
      let msg = canonical.iter().map(quote).join(", ");
      return make_error!("Canonical set contains 'gap' character: {msg}");
    }

    if canonical.contains(*unknown) {
      let msg = canonical.iter().map(quote).join(", ");
      return make_error!("Canonical set contains 'unknown' character: {msg}");
    }

    if ambiguous_keys.contains(*gap) {
      let msg = ambiguous.keys().map(quote).join(", ");
      return make_error!("Ambiguous set contains 'gap' character: {msg}");
    }

    if ambiguous_keys.contains(*gap) {
      let msg = ambiguous.keys().map(quote).join(", ");
      return make_error!("Ambiguous set contains 'unknown' character: {msg}");
    }

    {
      let ambig_gaps = ambiguous_set_map
        .iter()
        .map(|(key, vals)| (key, vals.difference(&canonical)))
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
  // use super::*;
  // use eyre::Report;
  // use indoc::indoc;
  // use pretty_assertions::assert_eq;

  // #[test]
  // fn test_alphabet_nuc() -> Result<(), Report> {
  //   let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
  //   let actual = json_write_str(&alphabet, JsonPretty(true))?;
  //   let expected = indoc! { /* language=json */ r#"{
  //     "all": [
  //       "-",
  //       "A",
  //       "B",
  //       "C",
  //       "D",
  //       "G",
  //       "H",
  //       "K",
  //       "M",
  //       "N",
  //       "R",
  //       "S",
  //       "T",
  //       "V",
  //       "W",
  //       "Y"
  //     ],
  //     "canonical": [
  //       "A",
  //       "C",
  //       "G",
  //       "T"
  //     ],
  //     "ambiguous": {
  //       "R": [
  //         "A",
  //         "G"
  //       ],
  //       "Y": [
  //         "C",
  //         "T"
  //       ],
  //       "S": [
  //         "C",
  //         "G"
  //       ],
  //       "W": [
  //         "A",
  //         "T"
  //       ],
  //       "K": [
  //         "G",
  //         "T"
  //       ],
  //       "M": [
  //         "A",
  //         "C"
  //       ],
  //       "D": [
  //         "A",
  //         "G",
  //         "T"
  //       ],
  //       "H": [
  //         "A",
  //         "C",
  //         "T"
  //       ],
  //       "B": [
  //         "C",
  //         "G",
  //         "T"
  //       ],
  //       "V": [
  //         "A",
  //         "C",
  //         "G"
  //       ]
  //     },
  //     "ambiguous_keys": [
  //       "B",
  //       "D",
  //       "H",
  //       "K",
  //       "M",
  //       "R",
  //       "S",
  //       "V",
  //       "W",
  //       "Y"
  //     ],
  //     "determined": [
  //       "A",
  //       "B",
  //       "C",
  //       "D",
  //       "G",
  //       "H",
  //       "K",
  //       "M",
  //       "R",
  //       "S",
  //       "T",
  //       "V",
  //       "W",
  //       "Y"
  //     ],
  //     "undetermined": [
  //       "-",
  //       "N"
  //     ],
  //     "unknown": "N",
  //     "gap": "-",
  //     "treat_gap_as_unknown": false,
  //     "profile_map": {
  //       "A": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "C": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "G": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0
  //         ]
  //       },
  //       "T": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0
  //         ]
  //       },
  //       "N": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0
  //         ]
  //       },
  //       "R": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           0.0,
  //           1.0,
  //           0.0
  //         ]
  //       },
  //       "Y": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           1.0,
  //           0.0,
  //           1.0
  //         ]
  //       },
  //       "S": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           1.0,
  //           1.0,
  //           0.0
  //         ]
  //       },
  //       "W": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           0.0,
  //           0.0,
  //           1.0
  //         ]
  //       },
  //       "K": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           1.0,
  //           1.0
  //         ]
  //       },
  //       "M": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           1.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "D": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           0.0,
  //           1.0,
  //           1.0
  //         ]
  //       },
  //       "H": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           1.0,
  //           0.0,
  //           1.0
  //         ]
  //       },
  //       "B": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           0.0,
  //           1.0,
  //           1.0,
  //           1.0
  //         ]
  //       },
  //       "V": {
  //         "v": 1,
  //         "dim": [
  //           4
  //         ],
  //         "data": [
  //           1.0,
  //           1.0,
  //           1.0,
  //           0.0
  //         ]
  //       }
  //     }
  //   }"#};
  //   assert_eq!(expected, actual);
  //   Ok(())
  // }
  //
  // #[test]
  // fn test_alphabet_aa() -> Result<(), Report> {
  //   let alphabet = Alphabet::new(AlphabetName::Aa, false)?;
  //   let actual = json_write_str(&alphabet, JsonPretty(true))?;
  //   let expected = indoc! { /* language=json */ r#"{
  //     "all": [
  //       "*",
  //       "-",
  //       "A",
  //       "B",
  //       "C",
  //       "D",
  //       "E",
  //       "F",
  //       "G",
  //       "H",
  //       "I",
  //       "J",
  //       "K",
  //       "L",
  //       "M",
  //       "N",
  //       "P",
  //       "Q",
  //       "R",
  //       "S",
  //       "T",
  //       "V",
  //       "W",
  //       "X",
  //       "Y",
  //       "Z"
  //     ],
  //     "canonical": [
  //       "*",
  //       "A",
  //       "C",
  //       "D",
  //       "E",
  //       "F",
  //       "G",
  //       "H",
  //       "I",
  //       "K",
  //       "L",
  //       "M",
  //       "N",
  //       "P",
  //       "Q",
  //       "R",
  //       "S",
  //       "T",
  //       "V",
  //       "W",
  //       "Y"
  //     ],
  //     "ambiguous": {
  //       "B": [
  //         "N",
  //         "D"
  //       ],
  //       "Z": [
  //         "Q",
  //         "E"
  //       ],
  //       "J": [
  //         "L",
  //         "I"
  //       ]
  //     },
  //     "ambiguous_keys": [
  //       "B",
  //       "J",
  //       "Z"
  //     ],
  //     "determined": [
  //       "*",
  //       "A",
  //       "B",
  //       "C",
  //       "D",
  //       "E",
  //       "F",
  //       "G",
  //       "H",
  //       "I",
  //       "J",
  //       "K",
  //       "L",
  //       "M",
  //       "N",
  //       "P",
  //       "Q",
  //       "R",
  //       "S",
  //       "T",
  //       "V",
  //       "W",
  //       "Y",
  //       "Z"
  //     ],
  //     "undetermined": [
  //       "-",
  //       "X"
  //     ],
  //     "unknown": "X",
  //     "gap": "-",
  //     "treat_gap_as_unknown": false,
  //     "profile_map": {
  //       "A": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "C": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "D": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "E": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "F": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "G": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "H": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "I": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "K": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "L": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "M": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "N": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "P": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "Q": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "R": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "S": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "T": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "V": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "W": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "Y": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0
  //         ]
  //       },
  //       "*": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0
  //         ]
  //       },
  //       "X": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0,
  //           1.0
  //         ]
  //       },
  //       "B": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "Z": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       },
  //       "J": {
  //         "v": 1,
  //         "dim": [
  //           21
  //         ],
  //         "data": [
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           1.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0,
  //           0.0
  //         ]
  //       }
  //     }
  //   }"#};
  //   assert_eq!(expected, actual);
  //   Ok(())
  // }
}
