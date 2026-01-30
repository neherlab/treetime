use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;

// Finds contiguous ranges (segments) in the sequence, such that for every character inside every range,
// the predicate function returns true and every range contains only the same letter.
//
// The predicate is a function that takes a character and returns boolean.
//
// For example if predicate returns `true` for characters A and C, this function will find ranges `AAAA`, `CCCCC`, `ACCCACAAA`
// but not `ZZZZZ`.
#[inline]
pub fn find_letter_ranges_by<F>(seq: &Seq, pred: F) -> Vec<(usize, usize)>
where
  F: Fn(AsciiChar) -> bool + Copy,
{
  let mut result = Vec::with_capacity(31);
  let mut i = 0;
  while i < seq.len() {
    if pred(seq[i]) {
      let start = i;
      i += 1;
      while i < seq.len() && pred(seq[i]) {
        i += 1;
      }
      result.push((start, i));
    } else {
      i += 1;
    }
  }
  result
}

pub mod old {
  #![allow(dead_code)]

  /// OLD, SLOWER VERSION
  pub fn find_letter_ranges_by(seq: &[char], pred: impl Fn(char) -> bool) -> Vec<(usize, usize)> {
    let len = seq.len();

    let mut result = vec![];
    let mut i = 0_usize;
    let mut start = 0_usize;
    let mut is_inside_range = false;
    while i < len {
      let letter = seq[i];

      // Find beginning of a range
      if pred(letter) {
        start = i;
        is_inside_range = true;
      }

      if is_inside_range {
        // Rewind forward until we find a mismatch
        while i < len && pred(seq[i]) {
          i += 1;
        }

        // We found the end of the current range, so now it's complete
        let end = i;

        // Remember the range
        result.push((start, end));

        is_inside_range = false;
      } else if i < len {
        i += 1;
      }
    }
    result
  }
}

/// Finds contiguous ranges (segments) consisting of a given letter in the sequence.
#[inline]
pub fn find_letter_ranges(seq: &Seq, letter: impl Into<AsciiChar> + Copy) -> Vec<(usize, usize)> {
  find_letter_ranges_by(seq, |candidate: AsciiChar| candidate == letter.into())
}

#[cfg(test)]
mod tests {
  use super::*;
  use rstest::rstest;

  #[rstest]
  #[case::empty_seq("",                 b'X', vec![])]
  #[case::no_match("GGATNNACA-ANTYGG", b'X', vec![])]
  #[case::all_match("XXX",              b'X', vec![(0, 3)])]
  #[case::middle_range("ATGXXXTTTT",       b'X', vec![(3, 6)])]
  #[case::two_ranges("ATGXXXCATGXXXXA",  b'X', vec![(3, 6), (10, 14)])]
  #[case::at_end("GCAXXXX",          b'X', vec![(3, 7)])]
  #[case::at_start("XXXXGCA",          b'X', vec![(0, 4)])]
  fn test_find_letter_ranges(#[case] seq: &str, #[case] letter: u8, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq.into(), letter);
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case::empty_seq("",                vec![])]
  #[case::no_ambiguous("GGATGCACATATTAGG",vec![])]
  #[case::all_ambiguous("NNN",             vec![(0, 3)])]
  #[case::middle_range("ATGNNNTTTT",      vec![(3, 6)])]
  #[case::two_ranges("ATGNNNCATGNNNNA", vec![(3, 6), (10, 14)])]
  #[case::at_end("GCANNNN",         vec![(3, 7)])]
  #[case::at_start("NNNNGCA",         vec![(0, 4)])]
  fn test_find_ambiguous_ranges(#[case] seq: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq.into(), b'N');
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case::empty_seq("",                vec![])]
  #[case::no_gaps("GGATGCACATATTAGG",vec![])]
  #[case::all_gaps("---",             vec![(0, 3)])]
  #[case::middle_range("ATG---TTTT",      vec![(3, 6)])]
  #[case::two_ranges("ATG---CATG----A", vec![(3, 6), (10, 14)])]
  #[case::at_end("GCA----",         vec![(3, 7)])]
  #[case::at_start("----GCA",         vec![(0, 4)])]
  fn test_find_gap_ranges(#[case] seq: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq.into(), b'-');
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case::empty_seq("",                 vec![])]
  #[case::no_undetermined("GGATGCACATATTAGG", vec![])]
  #[case::all_n("NNNN",             vec![(0, 4)])]
  #[case::all_gaps("----",             vec![(0, 4)])]
  #[case::n_then_gaps("NNNN----",         vec![(0, 8)])]
  #[case::gaps_then_n("----NNNN",         vec![(0, 8)])]
  #[case::gaps_n_gaps("----NNNN----",     vec![(0, 12)])]
  #[case::n_gaps_n("NNNN----NNNN",     vec![(0, 12)])]
  #[case::alternating_gaps_n("A----TNNNNG----C", vec![(1, 5), (6, 10), (11, 15)])]
  #[case::alternating_n_gaps("ANNNNT----GNNNNC", vec![(1, 5), (6, 10), (11, 15)])]
  #[case::n_middle_gaps_end("ATGNNNTTTT---",    vec![(3, 6), (10, 13)])]
  #[case::gaps_middle_n_end("ATG---TTTTNNN",    vec![(3, 6), (10, 13)])]
  fn test_find_undetermined_ranges(#[case] seq: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges_by(&seq.into(), |c| c == AsciiChar(b'N') || c == AsciiChar(b'-'));
    assert_eq!(expected, actual);
  }
}
