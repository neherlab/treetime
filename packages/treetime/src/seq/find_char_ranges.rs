use crate::seq::range_union::range_union;

// Finds contiguous ranges (segments) in the sequence, such that for every character inside every range,
// the predicate function returns true and every range contains only the same letter.
//
// The predicate is a function that takes a character and returns boolean.
//
// For example if predicate returns `true` for characters A and C, this function will find ranges `AAAA` and `CCCCC`,
// but not `ZZZ` or `ACCCAC`.
pub fn find_letter_ranges_by(seq: &[char], pred: impl Fn(char) -> bool) -> Vec<(usize, usize)> {
  let len = seq.len();

  let mut result = vec![];
  let mut i = 0_usize;
  let mut start = 0_usize;
  let mut found_maybe = Option::<char>::default();
  while i < len {
    let letter = seq[i];

    // Find beginning of a range
    if pred(letter) {
      start = i;
      found_maybe = Some(letter);
    }

    match found_maybe {
      // If there's a current range we are working on (for which we found a `start`), extend it
      Some(found) => {
        // Rewind forward until we find a mismatch
        while i < len && seq[i] == found {
          i += 1;
        }

        // We found the end of the current range, so now it's complete
        let end = i;

        // Remember the range
        result.push((start, end));

        found_maybe = None;
      }
      None => {
        if i < len {
          i += 1;
        }
      }
    }
  }
  result
}

/// Finds contiguous ranges (segments) consisting of a given letter in the sequence.
pub fn find_letter_ranges(seq: &[char], letter: char) -> Vec<(usize, usize)> {
  find_letter_ranges_by(seq, |candidate| candidate == letter)
}

pub fn find_ambiguous_ranges(seq: &[char]) -> Vec<(usize, usize)> {
  find_letter_ranges(seq, 'N')
}

pub fn find_gap_ranges(seq: &[char]) -> Vec<(usize, usize)> {
  find_letter_ranges(seq, '-')
}

pub fn find_undetermined_ranges(seq: &[char]) -> Vec<(usize, usize)> {
  range_union(&[find_letter_ranges_by(seq, |c| c == 'N' || c == '-')])
}

#[cfg(test)]
mod tests {
  use super::*;
  use rstest::rstest;

  fn to_char_array(seq: &str) -> Vec<char> {
    seq.chars().collect()
  }

  #[rstest]
  #[case("",                'X', vec![])]
  #[case("GGATNNACA-ANTYGG",'X', vec![])]
  #[case("XXX",             'X', vec![(0, 3)])]
  #[case("ATGXXXTTTT",      'X', vec![(3, 6)])]
  #[case("ATGXXXCATGXXXXA", 'X', vec![(3, 6), (10, 14)])]
  #[case("GCAXXXX",         'X', vec![(3, 7)])]
  #[case("XXXXGCA",         'X', vec![(0, 4)])]
  fn test_find_letter_ranges(#[case] seq: &str, #[case] letter: char, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&to_char_array(seq), letter);
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
    let actual = find_ambiguous_ranges(&to_char_array(seq));
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
    let actual = find_gap_ranges(&to_char_array(seq));
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
    let actual = find_undetermined_ranges(&to_char_array(seq));
    assert_eq!(expected, actual);
  }
}
