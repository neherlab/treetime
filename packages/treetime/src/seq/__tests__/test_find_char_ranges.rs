#[cfg(test)]
mod tests {
  use crate::seq::find_char_ranges::{find_letter_ranges, find_letter_ranges_by};
  use rstest::rstest;
  use treetime_primitives::{AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn seq(s: &str) -> Seq {
    Seq::try_from_str(s).unwrap()
  }

  #[rstest]
  #[case::empty_seq("",                 b'X', vec![])]
  #[case::no_match("GGATNNACA-ANTYGG", b'X', vec![])]
  #[case::all_match("XXX",              b'X', vec![(0, 3)])]
  #[case::middle_range("ATGXXXTTTT",       b'X', vec![(3, 6)])]
  #[case::two_ranges("ATGXXXCATGXXXXA",  b'X', vec![(3, 6), (10, 14)])]
  #[case::at_end("GCAXXXX",          b'X', vec![(3, 7)])]
  #[case::at_start("XXXXGCA",          b'X', vec![(0, 4)])]
  fn test_find_letter_ranges(#[case] s: &str, #[case] letter: u8, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq(s), c(letter));
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
  fn test_find_ambiguous_ranges(#[case] s: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq(s), c(b'N'));
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
  fn test_find_gap_ranges(#[case] s: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges(&seq(s), c(b'-'));
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
  fn test_find_undetermined_ranges(#[case] s: &str, #[case] expected: Vec<(usize, usize)>) {
    let actual = find_letter_ranges_by(&seq(s), |ch| ch == c(b'N') || ch == c(b'-'));
    assert_eq!(expected, actual);
  }
}
