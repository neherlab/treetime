use crate::make_error;
use eyre::Report;
use lazy_static::lazy_static;
use ndarray::{array, Array, Array1, Dimension};
use std::ops::Index;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alphabet {
  pub name: String,
  pub alphabet: Array1<char>,
  pub gap_index: Option<usize>,
  pub ambiguous: char,
}

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

impl Alphabet {
  pub fn new(name: &str) -> Result<Self, Report> {
    let alphabet = match name {
      "nuc" => Ok(ALPHABET_NUC.to_owned()),
      "nucleotide" => Ok(ALPHABET_NUC.to_owned()),
      "DNA" => Ok(ALPHABET_NUC.to_owned()),
      "nuc_nogap" => Ok(ALPHABET_NUC_NOGAP.to_owned()),
      "nucleotide_nogap" => Ok(ALPHABET_NUC_NOGAP.to_owned()),
      "DNA_nogap" => Ok(ALPHABET_NUC_NOGAP.to_owned()),
      "aa" => Ok(ALPHABET_AA.to_owned()),
      "aminoacid" => Ok(ALPHABET_AA.to_owned()),
      "aa_nogap" => Ok(ALPHABET_AA_NOGAP.to_owned()),
      "aminoacid_nogap" => Ok(ALPHABET_AA_NOGAP.to_owned()),
      _ => make_error!("Unknown alphabet: '{name}'"),
    }?;

    let ambiguous = match name {
      "nuc" => Ok(LETTER_AMBIGUOUS_NUC.to_owned()),
      "nucleotide" => Ok(LETTER_AMBIGUOUS_NUC.to_owned()),
      "DNA" => Ok(LETTER_AMBIGUOUS_NUC.to_owned()),
      "nuc_nogap" => Ok(LETTER_AMBIGUOUS_NUC_NOGAP.to_owned()),
      "nucleotide_nogap" => Ok(LETTER_AMBIGUOUS_NUC_NOGAP.to_owned()),
      "DNA_nogap" => Ok(LETTER_AMBIGUOUS_NUC_NOGAP.to_owned()),
      "aa" => Ok(LETTER_AMBIGUOUS_AA.to_owned()),
      "aminoacid" => Ok(LETTER_AMBIGUOUS_AA.to_owned()),
      "aa_nogap" => Ok(LETTER_AMBIGUOUS_AA_NOGAP.to_owned()),
      "aminoacid_nogap" => Ok(LETTER_AMBIGUOUS_AA_NOGAP.to_owned()),
      _ => make_error!("Unknown alphabet: '{name}'"),
    }?;

    let gap_index = alphabet.iter().position(|&x| x == '-');

    Ok(Self {
      name: name.to_owned(),
      alphabet,
      gap_index,
      ambiguous,
    })
  }

  pub fn indices_to_seq<D: Dimension>(&self, indices: &Array<usize, D>) -> Array<char, D> {
    indices.map(|&i| self.alphabet[i])
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
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn converts_indices_to_sequence() -> Result<(), Report> {
    let alphabet = Alphabet::new("nuc")?;
    assert_eq!(
      alphabet.indices_to_seq(&array![2, 3, 2, 4, 2, 2, 1]),
      array!['G', 'T', 'G', '-', 'G', 'G', 'C'],
    );
    Ok(())
  }

  #[rstest]
  fn converts_indices_to_sequence_2d() -> Result<(), Report> {
    let alphabet = Alphabet::new("nuc")?;
    assert_eq!(
      alphabet.indices_to_seq(&array![[2, 3, 2, 4, 2, 2, 1], [2, 3, 2, 4, 2, 2, 1]]),
      array![['G', 'T', 'G', '-', 'G', 'G', 'C'], ['G', 'T', 'G', '-', 'G', 'G', 'C']],
    );
    Ok(())
  }
}
