#[cfg(test)]
mod tests {
  use crate::fmt::string::{TruncateDirection, counted, name_list, truncate, truncate_right_with_ellipsis};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty_input(        ("",            5),  "")]
  #[case::shorter_than_limit( ("hello",       10), "hello")]
  #[case::exact_length(       ("hello world", 11), "hello world")]
  #[case::truncated(          ("hello world", 8),  "hello...")]
  #[case::only_ellipsis(      ("hello world", 3),  "...")]
  #[case::ellipsis_too_long(  ("hello",       2),  "he")]
  #[case::zero_limit(         ("hello",       0),  "")]
  #[trace]
  fn test_truncate_right_with_ellipsis_fits_limit(#[case] (input, max_len): (&str, usize), #[case] expected: &str) {
    let actual = truncate_right_with_ellipsis(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::fits(                     ("too short",   20, "...",    TruncateDirection::Right),  "too short")]
  #[case::right_empty_ellipsis(     ("hello world", 8,  "",       TruncateDirection::Right),  "hello wo")]
  #[case::right(                    ("hello world", 8,  "<cut>",  TruncateDirection::Right),  "hel<cut>")]
  #[case::left(                     ("foo,bar;baz", 8,  " >>",    TruncateDirection::Left),   " >>r;baz")]
  #[case::left_ellipsis_too_long(   ("xy",          1,  "...",    TruncateDirection::Left),   "y")]
  #[case::middle_even_remainder(    ("just over",   8,  "....",   TruncateDirection::Middle), "ju....er")]
  #[case::middle_odd_remainder(     ("hello world", 10, "[...]",  TruncateDirection::Middle), "he[...]rld")]
  #[case::middle_no_left_part(      ("abc",         2,  "z",      TruncateDirection::Middle), "zc")]
  #[case::middle_ellipsis_too_long( ("abcdef",      3,  "....",   TruncateDirection::Middle), "aef")]
  #[trace]
  fn test_truncate_with_custom_ellipsis_fits_limit(
    #[case] (input, max_len, ellipsis, direction): (&str, usize, &str, TruncateDirection),
    #[case] expected: &str,
  ) {
    let actual = truncate(input, max_len, Some(ellipsis), direction);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(       vec![],                                  "")]
  #[case::one(         vec!["A"],                               "A")]
  #[case::all_shown(   vec!["A", "B", "C"],                     "A, B, C")]
  #[case::one_more(    vec!["A", "B", "C", "D"],                "A, B, C, and 1 more")]
  #[case::many_more(   vec!["A", "B", "C", "D", "E", "F"],      "A, B, C, and 3 more")]
  #[trace]
  fn test_name_list_shows_first_three_and_counts_rest(#[case] names: Vec<&str>, #[case] expected: &str) {
    let names = names.into_iter().map(str::to_owned).collect::<Vec<_>>();
    assert_eq!(expected, name_list(&names));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero( 0, "0 tips")]
  #[case::one(  1, "1 tip")]
  #[case::many( 2, "2 tips")]
  #[trace]
  fn test_counted_pluralizes_noun_unless_one(#[case] count: usize, #[case] expected: &str) {
    assert_eq!(expected, counted(count, "tip"));
  }
}
