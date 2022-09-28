use clap::ArgEnum;
use eyre::Report;
use itertools::Itertools;
use lazy_static::lazy_static;
use ndarray::iter::Iter;
use ndarray::{array, Array, Array1, Dimension, Ix1};
use smart_default::SmartDefault;
use std::ops::Index;
use strum_macros::Display;

lazy_static! {
  static ref ALPHABET_NUC: Array1<char> = array!['A', 'C', 'G', 'T', '-'];
  static ref LETTER_AMBIGUOUS_NUC: char = 'N';
  static ref ALPHABET_NUC_NOGAP: Array1<char> = array!['A', 'C', 'G', 'T'];
  static ref LETTER_AMBIGUOUS_NUC_NOGAP: char = 'N';
  static ref ALPHABET_AA: Array1<char> = array![
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', '-'
  ];
  static ref LETTER_AMBIGUOUS_AA: char = 'X';
  static ref ALPHABET_AA_NOGAP: Array1<char> =
    array!['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
  static ref LETTER_AMBIGUOUS_AA_NOGAP: char = 'X';
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault, Display)]
#[clap(rename = "kebab-case")]
pub enum AlphabetName {
  #[default]
  Nuc,
  NucNogap,
  Aa,
  AaNogap,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alphabet {
  pub name: AlphabetName,
  pub alphabet: Array1<char>,
  pub gap_index: Option<usize>,
  pub ambiguous: char,
}

impl Alphabet {
  pub fn new(name: AlphabetName) -> Result<Self, Report> {
    let (alphabet, ambiguous) = match name {
      AlphabetName::Nuc => (ALPHABET_NUC.to_owned(), *LETTER_AMBIGUOUS_NUC),
      AlphabetName::NucNogap => (ALPHABET_NUC_NOGAP.to_owned(), *LETTER_AMBIGUOUS_NUC_NOGAP),
      AlphabetName::Aa => (ALPHABET_AA.to_owned(), *LETTER_AMBIGUOUS_AA),
      AlphabetName::AaNogap => (ALPHABET_AA_NOGAP.to_owned(), *LETTER_AMBIGUOUS_AA_NOGAP),
    };

    let gap_index = alphabet.iter().position(|&x| x == '-');

    Ok(Self {
      name,
      alphabet,
      gap_index,
      ambiguous,
    })
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
