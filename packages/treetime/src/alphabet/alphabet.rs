use clap::ArgEnum;
use eyre::Report;
use indexmap::IndexMap;
use itertools::Itertools;
use lazy_static::lazy_static;
use ndarray::iter::Iter;
use ndarray::{array, Array, Array1, Array2, Dimension, Ix1};
use smart_default::SmartDefault;
use strum_macros::Display;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault, Display)]
#[clap(rename = "kebab-case")]
pub enum AlphabetName {
  #[default]
  Nuc,
  NucNogap,
  Aa,
  AaNogap,
}

pub type ProfileMap = IndexMap<char, Array1<f64>>;

#[derive(Clone, Debug, PartialEq)]
pub struct Alphabet {
  pub alphabet: Array1<char>,
  pub profile_map: ProfileMap,
  pub gap_index: Option<usize>,
  pub ambiguous: char,
}

impl Alphabet {
  /// Creates one of the pre-defined alphabets
  pub fn new(name: AlphabetName) -> Result<Self, Report> {
    let (alphabet, profile_map, ambiguous) = match name {
      AlphabetName::Nuc => (
        ALPHABET_NUC.to_owned(),
        PROFILE_MAP_NUC.to_owned(),
        *LETTER_AMBIGUOUS_NUC,
      ),
      AlphabetName::NucNogap => (
        ALPHABET_NUC_NOGAP.to_owned(),
        PROFILE_MAP_NUC_NOGAP.to_owned(),
        *LETTER_AMBIGUOUS_NUC_NOGAP,
      ),
      AlphabetName::Aa => (
        ALPHABET_AA.to_owned(),
        PROFILE_MAP_NUC_NOGAP.to_owned(),
        *LETTER_AMBIGUOUS_AA,
      ),
      AlphabetName::AaNogap => (
        ALPHABET_AA_NOGAP.to_owned(),
        PROFILE_MAP_NUC_NOGAP.to_owned(),
        *LETTER_AMBIGUOUS_AA_NOGAP,
      ),
    };

    let gap_index = alphabet.iter().position(|&x| x == '-');

    Ok(Self {
      alphabet,
      profile_map,
      gap_index,
      ambiguous,
    })
  }

  /// Creates custom alphabet from a given set of letters and generates a trivial unambiguous profile map
  pub fn with_letters(letters: &[char], ambiguous: char) -> Result<Self, Report> {
    let eye = Array2::<f64>::eye(letters.len());

    let mut profile_map: ProfileMap = letters
      .iter()
      .zip(eye.rows())
      .map(|(s, x)| (*s, x.to_owned()))
      .collect();

    profile_map.insert(ambiguous, Array1::<f64>::ones(letters.len()));

    Ok(Self {
      alphabet: Array1::<char>::from_iter(letters.iter().copied()),
      profile_map,
      gap_index: None,
      ambiguous,
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

  pub fn indices_to_seq<D: Dimension>(&self, indices: &Array<usize, D>) -> Array<char, D> {
    indices.map(|&i| self.alphabet[i])
  }

  pub fn characters(&self) -> String {
    self.alphabet.iter().join("")
  }

  #[inline]
  pub fn len(&self) -> usize {
    self.alphabet.len()
  }

  #[inline]
  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  #[inline]
  pub const fn gap_index(&self) -> Option<usize> {
    self.gap_index
  }

  #[inline]
  pub const fn ambiguous(&self) -> char {
    self.ambiguous
  }

  #[allow(clippy::iter_without_into_iter)]
  pub fn iter(&self) -> Iter<'_, char, Ix1> {
    self.alphabet.iter()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn converts_indices_to_sequence() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    assert_eq!(
      alphabet.indices_to_seq(&array![2, 3, 2, 4, 2, 2, 1]),
      array!['G', 'T', 'G', '-', 'G', 'G', 'C'],
    );
    Ok(())
  }

  #[rstest]
  fn converts_indices_to_sequence_2d() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    assert_eq!(
      alphabet.indices_to_seq(&array![[2, 3, 2, 4, 2, 2, 1], [2, 3, 2, 4, 2, 2, 1]]),
      array![['G', 'T', 'G', '-', 'G', 'G', 'C'], ['G', 'T', 'G', '-', 'G', 'G', 'C']],
    );
    Ok(())
  }
}

lazy_static! {
  static ref ALPHABET_NUC: Array1<char> = array!['A', 'C', 'G', 'T', '-'];
  static ref LETTER_AMBIGUOUS_NUC: char = 'N';
  static ref PROFILE_MAP_NUC: ProfileMap = ProfileMap::from([
    ('A', array![1.0, 0.0, 0.0, 0.0, 0.0]),
    ('C', array![0.0, 1.0, 0.0, 0.0, 0.0]),
    ('G', array![0.0, 0.0, 1.0, 0.0, 0.0]),
    ('T', array![0.0, 0.0, 0.0, 1.0, 0.0]),
    ('-', array![0.0, 0.0, 0.0, 0.0, 1.0]),
    ('N', array![1.0, 1.0, 1.0, 1.0, 1.0]),
    ('X', array![1.0, 1.0, 1.0, 1.0, 1.0]),
    ('R', array![1.0, 0.0, 1.0, 0.0, 0.0]),
    ('Y', array![0.0, 1.0, 0.0, 1.0, 0.0]),
    ('S', array![0.0, 1.0, 1.0, 0.0, 0.0]),
    ('W', array![1.0, 0.0, 0.0, 1.0, 0.0]),
    ('K', array![0.0, 0.0, 1.0, 1.0, 0.0]),
    ('M', array![1.0, 1.0, 0.0, 0.0, 0.0]),
    ('D', array![1.0, 0.0, 1.0, 1.0, 0.0]),
    ('H', array![1.0, 1.0, 0.0, 1.0, 0.0]),
    ('B', array![0.0, 1.0, 1.0, 1.0, 0.0]),
    ('V', array![1.0, 1.0, 1.0, 0.0, 0.0]),
  ]);


  static ref ALPHABET_NUC_NOGAP: Array1<char> = array!['A', 'C', 'G', 'T'];
  static ref LETTER_AMBIGUOUS_NUC_NOGAP: char = 'N';
  static ref PROFILE_MAP_NUC_NOGAP: ProfileMap = ProfileMap::from([
    ('A', array![1.0, 0.0, 0.0, 0.0]),
    ('C', array![0.0, 1.0, 0.0, 0.0]),
    ('G', array![0.0, 0.0, 1.0, 0.0]),
    ('T', array![0.0, 0.0, 0.0, 1.0]),
    ('-', array![1.0, 1.0, 1.0, 1.0]), // gaps are completely ignored in distance computations
    ('N', array![1.0, 1.0, 1.0, 1.0]),
    ('X', array![1.0, 1.0, 1.0, 1.0]),
    ('R', array![1.0, 0.0, 1.0, 0.0]),
    ('Y', array![0.0, 1.0, 0.0, 1.0]),
    ('S', array![0.0, 1.0, 1.0, 0.0]),
    ('W', array![1.0, 0.0, 0.0, 1.0]),
    ('K', array![0.0, 0.0, 1.0, 1.0]),
    ('M', array![1.0, 1.0, 0.0, 0.0]),
    ('D', array![1.0, 0.0, 1.0, 1.0]),
    ('H', array![1.0, 1.0, 0.0, 1.0]),
    ('B', array![0.0, 1.0, 1.0, 1.0]),
    ('V', array![1.0, 1.0, 1.0, 0.0]),
  ]);


  static ref ALPHABET_AA: Array1<char> = array![
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', '-'
  ];
  static ref LETTER_AMBIGUOUS_AA: char = 'X';
  static ref PROFILE_MAP_AA: ProfileMap = ProfileMap::from([
    ('A', array![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Alanine         Ala
    ('C', array![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Cysteine        Cys
    ('D', array![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Aspartic AciD   Asp
    ('E', array![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glutamic Acid   Glu
    ('F', array![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Phenylalanine   Phe
    ('G', array![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glycine         Gly
    ('H', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Histidine       His
    ('I', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Isoleucine      Ile
    ('K', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Lysine          Lys
    ('L', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Leucine         Leu
    ('M', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Methionine      Met
    ('N', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // AsparagiNe      Asn
    ('P', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Proline         Pro
    ('Q', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glutamine       Gln
    ('R', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // ARginine        Arg
    ('S', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Serine          Ser
    ('T', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Threonine       Thr
    ('V', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]), // Valine          Val
    ('W', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]), // Tryptophan      Trp
    ('Y', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]), // Tyrosine        Tyr
    ('*', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]), // stop
    ('-', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]), // gap
    ('X', array![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]), // not specified/any
    ('B', array![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Asparagine/Aspartic Acid    Asx
    ('Z', array![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glutamine/Glutamic Acid     Glx
  ]);


  static ref ALPHABET_AA_NOGAP: Array1<char> =
    array!['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
  static ref LETTER_AMBIGUOUS_AA_NOGAP: char = 'X';
  static ref PROFILE_MAP_AA_NOGAP: ProfileMap = ProfileMap::from([
     ('A', array![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Alanine         Ala
     ('C', array![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Cysteine        Cys
     ('D', array![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Aspartic AciD   Asp
     ('E', array![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glutamic Acid   Glu
     ('F', array![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Phenylalanine   Phe
     ('G', array![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glycine         Gly
     ('H', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Histidine       His
     ('I', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Isoleucine      Ile
     ('K', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Lysine          Lys
     ('L', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Leucine         Leu
     ('M', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Methionine      Met
     ('N', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // AsparagiNe      Asn
     ('P', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Proline         Pro
     ('Q', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glutamine       Gln
     ('R', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // ARginine        Arg
     ('S', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]), // Serine          Ser
     ('T', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]), // Threonine       Thr
     ('V', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]), // Valine          Val
     ('W', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]), // Tryptophan      Trp
     ('Y', array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]), // Tyrosine        Tyr
     ('X', array![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]), // not specified/any
     ('B', array![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Asparagine/Aspartic Acid    Asx
     ('Z', array![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), // Glutamine/Glutamic Acid     Glx
  ]);
}
