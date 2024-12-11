use crate::representation::seq::Seq;

// Finds contiguous ranges (segments) in the sequence, such that for every character inside every range,
// the predicate function returns true and every range contains only the same letter.
//
// The predicate is a function that takes a character and returns boolean.
//
// For example if predicate returns `true` for characters A and C, this function will find ranges `AAAA`, `CCCCC`, `ACCCACAAA`
// but not `ZZZZZ`.
#[inline]
pub fn find_letter_ranges_by(seq: &Seq, pred: impl Fn(u8) -> bool + Copy) -> Vec<(usize, usize)> {
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
pub fn find_letter_ranges(seq: &Seq, letter: u8) -> Vec<(usize, usize)> {
  find_letter_ranges_by(seq, |candidate| candidate == letter)
}

#[cfg(test)]
mod tests {
  use super::*;
  use rstest::rstest;

  #[rstest]
  #[case("",                 b'X', vec![])]
  #[case("GGATNNACA-ANTYGG", b'X', vec![])]
  #[case("XXX",              b'X', vec![(0, 3)])]
  #[case("ATGXXXTTTT",       b'X', vec![(3, 6)])]
  #[case("ATGXXXCATGXXXXA",  b'X', vec![(3, 6), (10, 14)])]
  #[case("GCAXXXX",          b'X', vec![(3, 7)])]
  #[case("XXXXGCA",          b'X', vec![(0, 4)])]
  fn test_find_letter_ranges(#[case] seq: &str, #[case] letter: u8, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq.into(), letter);
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case("",                vec![])]
  #[case("GGATGCACATATTAGG",vec![])]
  #[case("NNN",             vec![(0, 3)])]
  #[case("ATGNNNTTTT",      vec![(3, 6)])]
  #[case("ATGNNNCATGNNNNA", vec![(3, 6), (10, 14)])]
  #[case("GCANNNN",         vec![(3, 7)])]
  #[case("NNNNGCA",         vec![(0, 4)])]
  fn test_find_ambiguous_ranges(#[case] seq: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq.into(), b'N');
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case("",                vec![])]
  #[case("GGATGCACATATTAGG",vec![])]
  #[case("---",             vec![(0, 3)])]
  #[case("ATG---TTTT",      vec![(3, 6)])]
  #[case("ATG---CATG----A", vec![(3, 6), (10, 14)])]
  #[case("GCA----",         vec![(3, 7)])]
  #[case("----GCA",         vec![(0, 4)])]
  fn test_find_gap_ranges(#[case] seq: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq.into(), b'-');
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case("",                 vec![])]
  #[case("GGATGCACATATTAGG", vec![])]
  #[case("NNNN",             vec![(0, 4)])]
  #[case("----",             vec![(0, 4)])]
  #[case("NNNN----",         vec![(0, 8)])]
  #[case("----NNNN",         vec![(0, 8)])]
  #[case("----NNNN----",     vec![(0, 12)])]
  #[case("NNNN----NNNN",     vec![(0, 12)])]
  #[case("A----TNNNNG----C", vec![(1, 5), (6, 10), (11, 15)])]
  #[case("ANNNNT----GNNNNC", vec![(1, 5), (6, 10), (11, 15)])]
  #[case("ATGNNNTTTT---",    vec![(3, 6), (10, 13)])]
  #[case("ATG---TTTTNNN",    vec![(3, 6), (10, 13)])]
  fn test_find_undetermined_ranges(#[case] seq: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges_by(&seq.into(), |c| c == b'N' || c == b'-');
    assert_eq!(expected, actual);
  }
}
