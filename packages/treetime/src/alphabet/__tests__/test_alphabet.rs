#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName, FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
  use crate::alphabet::alphabet_config::AlphabetConfig;
  use crate::pretty_assert_ulps_eq;
  use crate::vec_u8;
  use indexmap::indexmap;
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::{AlphabetLike, AsciiChar};

  #[test]
  fn test_alphabet_default() {
    let alphabet = Alphabet::default();
    let expected_n_canonical = 4;
    assert_eq!(expected_n_canonical, alphabet.n_canonical());
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc, 4)]
  #[case::aa(AlphabetName::Aa, 21)]
  #[trace]
  fn test_alphabet_new_n_canonical(#[case] name: AlphabetName, #[case] expected: usize) {
    let alphabet = Alphabet::new(name).unwrap();
    assert_eq!(expected, alphabet.n_canonical());
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc, 10)]
  #[case::aa(AlphabetName::Aa, 3)]
  #[trace]
  fn test_alphabet_new_n_ambiguous(#[case] name: AlphabetName, #[case] expected: usize) {
    let alphabet = Alphabet::new(name).unwrap();
    assert_eq!(expected, alphabet.n_ambiguous());
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc, 2)]
  #[case::aa(AlphabetName::Aa, 2)]
  #[trace]
  fn test_alphabet_new_n_undetermined(#[case] name: AlphabetName, #[case] expected: usize) {
    let alphabet = Alphabet::new(name).unwrap();
    assert_eq!(expected, alphabet.n_undetermined());
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc, b'N')]
  #[case::aa(AlphabetName::Aa, b'X')]
  #[trace]
  fn test_alphabet_unknown(#[case] name: AlphabetName, #[case] expected: u8) {
    let alphabet = Alphabet::new(name).unwrap();
    assert_eq!(AsciiChar::from_byte_unchecked(expected), alphabet.unknown());
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc)]
  #[case::aa(AlphabetName::Aa)]
  #[trace]
  fn test_alphabet_gap(#[case] name: AlphabetName) {
    let alphabet = Alphabet::new(name).unwrap();
    let expected = AsciiChar::from_byte_unchecked(b'-');
    assert_eq!(expected, alphabet.gap());
  }

  #[test]
  fn test_alphabet_nuc_canonical_chars() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let actual = alphabet.canonical().collect_vec();
    let expected = vec![
      AsciiChar::from_byte_unchecked(b'A'),
      AsciiChar::from_byte_unchecked(b'C'),
      AsciiChar::from_byte_unchecked(b'G'),
      AsciiChar::from_byte_unchecked(b'T'),
    ];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_alphabet_aa_canonical_chars() {
    let alphabet = Alphabet::new(AlphabetName::Aa).unwrap();
    let actual = alphabet.canonical().collect_vec();
    let expected = vec![
      AsciiChar::from_byte_unchecked(b'*'),
      AsciiChar::from_byte_unchecked(b'A'),
      AsciiChar::from_byte_unchecked(b'C'),
      AsciiChar::from_byte_unchecked(b'D'),
      AsciiChar::from_byte_unchecked(b'E'),
      AsciiChar::from_byte_unchecked(b'F'),
      AsciiChar::from_byte_unchecked(b'G'),
      AsciiChar::from_byte_unchecked(b'H'),
      AsciiChar::from_byte_unchecked(b'I'),
      AsciiChar::from_byte_unchecked(b'K'),
      AsciiChar::from_byte_unchecked(b'L'),
      AsciiChar::from_byte_unchecked(b'M'),
      AsciiChar::from_byte_unchecked(b'N'),
      AsciiChar::from_byte_unchecked(b'P'),
      AsciiChar::from_byte_unchecked(b'Q'),
      AsciiChar::from_byte_unchecked(b'R'),
      AsciiChar::from_byte_unchecked(b'S'),
      AsciiChar::from_byte_unchecked(b'T'),
      AsciiChar::from_byte_unchecked(b'V'),
      AsciiChar::from_byte_unchecked(b'W'),
      AsciiChar::from_byte_unchecked(b'Y'),
    ];
    assert_eq!(expected, actual);
  }

  #[rstest]
  #[case::a(b'A', true)]
  #[case::c(b'C', true)]
  #[case::g(b'G', true)]
  #[case::t(b'T', true)]
  #[case::n(b'N', false)]
  #[case::gap(b'-', false)]
  #[case::r(b'R', false)]
  #[trace]
  fn test_alphabet_nuc_is_canonical(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.is_canonical(AsciiChar::from_byte_unchecked(c)));
  }

  #[rstest]
  #[case::r(b'R', true)]
  #[case::y(b'Y', true)]
  #[case::s(b'S', true)]
  #[case::w(b'W', true)]
  #[case::a(b'A', false)]
  #[case::n(b'N', false)]
  #[case::gap(b'-', false)]
  #[trace]
  fn test_alphabet_nuc_is_ambiguous(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.is_ambiguous(AsciiChar::from_byte_unchecked(c)));
  }

  #[rstest]
  #[case::a(b'A', true)]
  #[case::r(b'R', true)]
  #[case::n(b'N', false)]
  #[case::gap(b'-', false)]
  #[trace]
  fn test_alphabet_nuc_is_determined(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.is_determined(AsciiChar::from_byte_unchecked(c)));
  }

  #[rstest]
  #[case::n(b'N', true)]
  #[case::gap(b'-', true)]
  #[case::a(b'A', false)]
  #[case::r(b'R', false)]
  #[trace]
  fn test_alphabet_nuc_is_undetermined(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.is_undetermined(AsciiChar::from_byte_unchecked(c)));
  }

  #[rstest]
  #[case::n(b'N', true)]
  #[case::gap(b'-', false)]
  #[case::a(b'A', false)]
  #[trace]
  fn test_alphabet_nuc_is_unknown(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.is_unknown(AsciiChar::from_byte_unchecked(c)));
  }

  #[rstest]
  #[case::gap(b'-', true)]
  #[case::n(b'N', false)]
  #[case::a(b'A', false)]
  #[trace]
  fn test_alphabet_nuc_is_gap(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.is_gap(AsciiChar::from_byte_unchecked(c)));
  }

  #[rstest]
  #[case::a(b'A', true)]
  #[case::n(b'N', true)]
  #[case::gap(b'-', true)]
  #[case::r(b'R', true)]
  #[case::x(b'X', false)]
  #[case::z(b'Z', false)]
  #[trace]
  fn test_alphabet_nuc_contains(#[case] c: u8, #[case] expected: bool) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.contains(AsciiChar::from_byte_unchecked(c)));
  }

  #[test]
  fn test_alphabet_nuc_n_chars() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let expected = 16;
    assert_eq!(expected, alphabet.n_chars());
  }

  #[test]
  fn test_alphabet_nuc_n_determined() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let expected = 14;
    assert_eq!(expected, alphabet.n_determined());
  }

  #[rstest]
  #[case::a(b'A', 0)]
  #[case::c(b'C', 1)]
  #[case::g(b'G', 2)]
  #[case::t(b'T', 3)]
  #[trace]
  fn test_alphabet_nuc_index(#[case] c: u8, #[case] expected: usize) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(expected, alphabet.index(c as usize).unwrap());
  }

  #[rstest]
  #[case::idx0(0, b'A')]
  #[case::idx1(1, b'C')]
  #[case::idx2(2, b'G')]
  #[case::idx3(3, b'T')]
  #[trace]
  fn test_alphabet_nuc_char(#[case] index: usize, #[case] expected: u8) {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert_eq!(AsciiChar::from_byte_unchecked(expected), alphabet.char(index));
  }

  #[test]
  fn test_alphabet_nuc_get_profile_canonical() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

    let profile_a = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'A')).unwrap();
    let expected_a = array![1.0, 0.0, 0.0, 0.0];
    pretty_assert_ulps_eq!(expected_a, profile_a, max_ulps = 4);

    let profile_c = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'C')).unwrap();
    let expected_c = array![0.0, 1.0, 0.0, 0.0];
    pretty_assert_ulps_eq!(expected_c, profile_c, max_ulps = 4);

    let profile_g = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'G')).unwrap();
    let expected_g = array![0.0, 0.0, 1.0, 0.0];
    pretty_assert_ulps_eq!(expected_g, profile_g, max_ulps = 4);

    let profile_t = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'T')).unwrap();
    let expected_t = array![0.0, 0.0, 0.0, 1.0];
    pretty_assert_ulps_eq!(expected_t, profile_t, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_nuc_get_profile_ambiguous() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

    let profile_r = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'R')).unwrap();
    let expected_r = array![1.0, 0.0, 1.0, 0.0];
    pretty_assert_ulps_eq!(expected_r, profile_r, max_ulps = 4);

    let profile_y = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'Y')).unwrap();
    let expected_y = array![0.0, 1.0, 0.0, 1.0];
    pretty_assert_ulps_eq!(expected_y, profile_y, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_nuc_get_profile_unknown() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let profile_n = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'N')).unwrap();
    let expected_n = array![1.0, 1.0, 1.0, 1.0];
    pretty_assert_ulps_eq!(expected_n, profile_n, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_nuc_get_code() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

    let actual = alphabet.get_code(&array![1.0, 0.0, 0.0, 0.0]).unwrap();
    let expected = AsciiChar::from_byte_unchecked(b'A');
    assert_eq!(expected, actual);

    let actual = alphabet.get_code(&array![0.0, 1.0, 0.0, 1.0]).unwrap();
    let expected = AsciiChar::from_byte_unchecked(b'Y');
    assert_eq!(expected, actual);

    let actual = alphabet.get_code(&array![1.0, 1.0, 1.0, 1.0]).unwrap();
    let expected = AsciiChar::from_byte_unchecked(b'N');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_alphabet_nuc_char_to_set_canonical() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let set = alphabet.char_to_set(AsciiChar::from_byte_unchecked(b'A'));
    let expected_len = 1;
    assert_eq!(expected_len, set.len());
    assert!(set.contains(AsciiChar::from_byte_unchecked(b'A')));
  }

  #[test]
  fn test_alphabet_nuc_char_to_set_ambiguous() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let set = alphabet.char_to_set(AsciiChar::from_byte_unchecked(b'R'));
    let expected_len = 2;
    assert_eq!(expected_len, set.len());
    assert!(set.contains(AsciiChar::from_byte_unchecked(b'A')));
    assert!(set.contains(AsciiChar::from_byte_unchecked(b'G')));
  }

  #[test]
  fn test_alphabet_nuc_set_to_char() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let set = alphabet.char_to_set(AsciiChar::from_byte_unchecked(b'R'));
    let actual = alphabet.set_to_char(set);
    let expected = AsciiChar::from_byte_unchecked(b'R');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_alphabet_nuc_construct_profile() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let actual = alphabet
      .construct_profile([
        AsciiChar::from_byte_unchecked(b'A'),
        AsciiChar::from_byte_unchecked(b'G'),
      ])
      .unwrap();
    let expected = array![1.0, 0.0, 1.0, 0.0];
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_nuc_seq2prof() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let seq = [
      AsciiChar::from_byte_unchecked(b'A'),
      AsciiChar::from_byte_unchecked(b'C'),
      AsciiChar::from_byte_unchecked(b'G'),
    ];
    let actual = alphabet.seq2prof(&seq).unwrap();
    let expected_shape = [3, 4];
    assert_eq!(expected_shape, actual.shape());

    let expected_row0 = array![1.0, 0.0, 0.0, 0.0];
    pretty_assert_ulps_eq!(expected_row0, actual.row(0).to_owned(), max_ulps = 4);
  }

  #[test]
  fn test_alphabet_gap_is_not_unknown() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    assert!(!alphabet.is_unknown(AsciiChar::from_byte_unchecked(b'-')));
    assert!(alphabet.is_gap(AsciiChar::from_byte_unchecked(b'-')));
  }

  #[test]
  fn test_alphabet_gap_profile_matches_unknown() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let profile_gap = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'-')).unwrap();
    let profile_unknown = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'N')).unwrap();
    pretty_assert_ulps_eq!(profile_unknown, profile_gap, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_reserved_constants() {
    let expected_non_char = AsciiChar::from_byte_unchecked(b'.');
    let expected_variable_char = AsciiChar::from_byte_unchecked(b'~');
    let expected_fill_char = AsciiChar::from_byte_unchecked(b' ');
    assert_eq!(expected_non_char, NON_CHAR);
    assert_eq!(expected_variable_char, VARIABLE_CHAR);
    assert_eq!(expected_fill_char, FILL_CHAR);
  }

  #[test]
  fn test_alphabet_aa_ambiguous() {
    let alphabet = Alphabet::new(AlphabetName::Aa).unwrap();
    let ambiguous = alphabet.ambiguous().collect_vec();
    let expected = vec![
      AsciiChar::from_byte_unchecked(b'B'),
      AsciiChar::from_byte_unchecked(b'J'),
      AsciiChar::from_byte_unchecked(b'Z'),
    ];
    assert_eq!(expected, ambiguous);
  }

  #[test]
  fn test_alphabet_aa_b_maps_to_nd() {
    let alphabet = Alphabet::new(AlphabetName::Aa).unwrap();
    let set = alphabet.char_to_set(AsciiChar::from_byte_unchecked(b'B'));
    assert!(set.contains(AsciiChar::from_byte_unchecked(b'N')));
    assert!(set.contains(AsciiChar::from_byte_unchecked(b'D')));
    let expected_len = 2;
    assert_eq!(expected_len, set.len());
  }

  #[test]
  fn test_alphabet_with_config_empty_canonical() {
    let config = AlphabetConfig {
      canonical: vec![],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    let result = Alphabet::with_config(&config);
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("empty"));
  }

  #[test]
  fn test_alphabet_nuc_chars_iterator() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let chars = alphabet.chars().collect_vec();
    let expected_len = 16;
    assert_eq!(expected_len, chars.len());
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'A')));
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'C')));
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'G')));
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'T')));
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'N')));
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'-')));
    assert!(chars.contains(&AsciiChar::from_byte_unchecked(b'R')));
  }

  #[test]
  fn test_alphabet_nuc_determined_iterator() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let determined = alphabet.determined().collect_vec();
    let expected_len = 14;
    assert_eq!(expected_len, determined.len());
    assert!(determined.contains(&AsciiChar::from_byte_unchecked(b'A')));
    assert!(determined.contains(&AsciiChar::from_byte_unchecked(b'R')));
    assert!(!determined.contains(&AsciiChar::from_byte_unchecked(b'N')));
    assert!(!determined.contains(&AsciiChar::from_byte_unchecked(b'-')));
  }

  #[test]
  fn test_alphabet_nuc_undetermined_iterator() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let undetermined = alphabet.undetermined().collect_vec();
    let expected_len = 2;
    assert_eq!(expected_len, undetermined.len());
    assert!(undetermined.contains(&AsciiChar::from_byte_unchecked(b'N')));
    assert!(undetermined.contains(&AsciiChar::from_byte_unchecked(b'-')));
  }

  #[test]
  fn test_alphabet_nuc_ambiguous_iterator() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let ambiguous = alphabet.ambiguous().collect_vec();
    let expected_len = 10;
    assert_eq!(expected_len, ambiguous.len());
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'R')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'Y')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'S')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'W')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'K')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'M')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'D')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'H')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'B')));
    assert!(ambiguous.contains(&AsciiChar::from_byte_unchecked(b'V')));
  }

  #[test]
  fn test_alphabet_construct_profile_with_ambiguous() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let actual = alphabet
      .construct_profile([AsciiChar::from_byte_unchecked(b'R')])
      .unwrap();
    let expected = array![1.0, 0.0, 1.0, 0.0];
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_construct_profile_empty() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let actual = alphabet.construct_profile::<[AsciiChar; 0], _>([]).unwrap();
    let expected = array![0.0, 0.0, 0.0, 0.0];
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
  }

  #[test]
  fn test_alphabet_seq2prof_with_ambiguous() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let seq = [
      AsciiChar::from_byte_unchecked(b'R'),
      AsciiChar::from_byte_unchecked(b'Y'),
    ];
    let actual = alphabet.seq2prof(&seq).unwrap();
    let expected_shape = [2, 4];
    assert_eq!(expected_shape, actual.shape());

    let expected_row0 = array![1.0, 0.0, 1.0, 0.0];
    pretty_assert_ulps_eq!(expected_row0, actual.row(0).to_owned(), max_ulps = 4);

    let expected_row1 = array![0.0, 1.0, 0.0, 1.0];
    pretty_assert_ulps_eq!(expected_row1, actual.row(1).to_owned(), max_ulps = 4);
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc, "Nuc")]
  #[case::aa(AlphabetName::Aa, "Aa")]
  #[trace]
  fn test_alphabet_name_display(#[case] name: AlphabetName, #[case] expected: &str) {
    let actual = name.to_string();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_alphabet_name_default() {
    let expected = AlphabetName::Nuc;
    let actual = AlphabetName::default();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_alphabet_with_custom_config() {
    let config = AlphabetConfig {
      canonical: vec_u8!['X', 'Y', 'Z'],
      ambiguous: indexmap! {
        b'W' => vec_u8!['X', 'Y'],
      },
      unknown: b'?',
      gap: b'-',
    };
    let alphabet = Alphabet::with_config(&config).unwrap();

    let expected_n_canonical = 3;
    assert_eq!(expected_n_canonical, alphabet.n_canonical());

    let expected_n_ambiguous = 1;
    assert_eq!(expected_n_ambiguous, alphabet.n_ambiguous());

    let expected_unknown = AsciiChar::from_byte_unchecked(b'?');
    assert_eq!(expected_unknown, alphabet.unknown());

    let expected_gap = AsciiChar::from_byte_unchecked(b'-');
    assert_eq!(expected_gap, alphabet.gap());

    let profile_w = alphabet.get_profile(AsciiChar::from_byte_unchecked(b'W')).unwrap();
    let expected_w = array![1.0, 1.0, 0.0];
    pretty_assert_ulps_eq!(expected_w, profile_w, max_ulps = 4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::r(AlphabetName::Nuc, b'R', vec![b'A', b'G'])]
  #[case::y(AlphabetName::Nuc, b'Y', vec![b'C', b'T'])]
  #[case::s(AlphabetName::Nuc, b'S', vec![b'C', b'G'])]
  #[case::w(AlphabetName::Nuc, b'W', vec![b'A', b'T'])]
  #[case::k(AlphabetName::Nuc, b'K', vec![b'G', b'T'])]
  #[case::m(AlphabetName::Nuc, b'M', vec![b'A', b'C'])]
  #[case::d(AlphabetName::Nuc, b'D', vec![b'A', b'G', b'T'])]
  #[case::h(AlphabetName::Nuc, b'H', vec![b'A', b'C', b'T'])]
  #[case::b(AlphabetName::Nuc, b'B', vec![b'C', b'G', b'T'])]
  #[case::v(AlphabetName::Nuc, b'V', vec![b'A', b'C', b'G'])]
  #[case::aa_b(AlphabetName::Aa, b'B', vec![b'N', b'D'])]
  #[case::aa_z(AlphabetName::Aa, b'Z', vec![b'Q', b'E'])]
  #[case::aa_j(AlphabetName::Aa, b'J', vec![b'L', b'I'])]
  #[trace]
  fn test_alphabet_ambiguous_mapping(
    #[case] name: AlphabetName,
    #[case] ambig: u8,
    #[case] expected_chars: Vec<u8>,
  ) {
    let alphabet = Alphabet::new(name).unwrap();
    let set = alphabet.char_to_set(AsciiChar::from_byte_unchecked(ambig));
    assert_eq!(expected_chars.len(), set.len());
    for c in expected_chars {
      assert!(set.contains(AsciiChar::from_byte_unchecked(c)));
    }
  }

  #[test]
  fn test_alphabet_char_index_roundtrip() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    for i in 0..alphabet.n_canonical() {
      let c = alphabet.char(i);
      let idx = alphabet.index(c).unwrap();
      assert_eq!(i, idx);
    }
  }

  #[test]
  fn test_alphabet_set_to_char_canonical_roundtrip() {
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    for c in alphabet.canonical() {
      let set = alphabet.char_to_set(c);
      let back = alphabet.set_to_char(set);
      assert_eq!(c, back);
    }
  }

  #[rstest]
  #[case::nuc(AlphabetName::Nuc)]
  #[case::aa(AlphabetName::Aa)]
  #[case::aa_no_stop(AlphabetName::AaNoStop)]
  #[trace]
  fn test_alphabet_profile_index_consistency(#[case] name: AlphabetName) {
    let alphabet = Alphabet::new(name).unwrap();
    for c in alphabet.canonical() {
      let index = alphabet.index(c).unwrap();
      let profile = alphabet.get_profile(c).unwrap();
      pretty_assert_ulps_eq!(1.0, profile[index], max_ulps = 4);

      for j in 0..alphabet.n_canonical() {
        if j != index {
          pretty_assert_ulps_eq!(0.0, profile[j], max_ulps = 4);
        }
      }

      let roundtrip = alphabet.char(index);
      assert_eq!(c, roundtrip, "char(index({c})) should roundtrip");
    }
  }

  #[test]
  fn test_alphabet_aa_star_index_matches_profile() {
    let alphabet = Alphabet::new(AlphabetName::Aa).unwrap();
    let star = AsciiChar::from_byte_unchecked(b'*');
    let index = alphabet.index(star).unwrap();
    let profile = alphabet.get_profile(star).unwrap();
    pretty_assert_ulps_eq!(1.0, profile[index], max_ulps = 4);
  }

  #[test]
  fn test_alphabet_char_index_roundtrip_aa() {
    let alphabet = Alphabet::new(AlphabetName::Aa).unwrap();
    for i in 0..alphabet.n_canonical() {
      let c = alphabet.char(i);
      let idx = alphabet.index(c).unwrap();
      assert_eq!(i, idx);
    }
  }
}
