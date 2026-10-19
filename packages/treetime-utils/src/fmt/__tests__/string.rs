#[cfg(test)]
mod tests {
  use crate::fmt::string::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::basic_truncate("hello world",                     8, "hello...")]
  #[case::short_limit("hello world",                     4, "h...")]
  #[case::exact_length("hello world",                    11, "hello world")]
  #[case::no_truncate_needed("hello",                          10, "hello")]
  #[case::five_chars("hello world",                     5, "he...")]
  #[case::only_ellipsis("hello world",                     3, "...")]
  #[case::empty_input("",                                5, "")]
  #[case::partial("abcdefghij",                      7, "abcd...")]
  #[case::zero_limit("hello",                           0, "")]
  #[case::two_chars_no_ellipsis("hello",                           2, "he")]
  #[case::one_char_no_ellipsis("hello",                           1, "h")]
  #[case::exact_fit("hello",                           5, "hello")]
  #[case::longer_limit("hello",                           6, "hello")]
  #[case::special_chars_semicolon("foo,bar;baz:qux",                10, "foo,bar...")]
  #[case::email_address("test@example.com",               12, "test@exam...")]
  #[case::file_path("path/to/file.txt",               12, "path/to/f...")]
  #[case::has_ellipsis("...already...has...ellipsis...", 20, "...already...has....")]
  #[case::symbols("!@#$%^&*()",                      7, "!@#$...")]
  #[case::with_spaces("a b c d e f g",                  10, "a b c d...")]
  #[trace]
  fn test_truncate_right_with_ellipsis(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] expected: &str,
  ) {
    let actual = truncate_right_with_ellipsis(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right_empty_ellipsis("hello world",          8, "",       TruncateDirection::Right,  "hello wo")]
  #[case::right_custom_cut("hello world",          8, "<cut>",  TruncateDirection::Right,  "hel<cut>")]
  #[case::middle_brackets("hello world",         10, "[...]",  TruncateDirection::Middle, "he[...]rld")]
  #[case::left_arrow("foo,bar;baz",          8, " >>",    TruncateDirection::Left,   " >>r;baz")]
  #[case::left_double_arrow("test@example.com",    12, ">>",     TruncateDirection::Left,   ">>xample.com")]
  #[case::middle_snip("path/to/file.txt",    15, "<snip>", TruncateDirection::Middle, "path<snip>e.txt")]
  #[case::middle_tilde("a b c d e f g",       10, " ~ ",    TruncateDirection::Middle, "a b ~  f g")]
  #[case::right_asterisks("!@#$%^&*()",           7, "***",    TruncateDirection::Right,  "!@#$***")]
  #[case::no_truncate_needed("too short",           20, "...",    TruncateDirection::Right,  "too short")]
  #[case::exact_length("exact",                5, "...",    TruncateDirection::Right,  "exact")]
  #[case::middle_four_dots("just over",            8, "....",   TruncateDirection::Middle, "ju....er")]
  #[case::input_has_ellipsis("ellipsis in input...", 15, "...",   TruncateDirection::Right,  "ellipsis in ...")]
  #[case::middle_angle_brackets("multiple words here", 12, " <> ",   TruncateDirection::Middle, "mult <> here")]
  #[case::single_char_input("x",                    1, "...",    TruncateDirection::Right,  "x")]
  #[case::two_char_input_left("xy",                   1, "...",    TruncateDirection::Left,   "y")]
  #[case::single_char_ellipsis("abc",                  2, "z",      TruncateDirection::Middle, "zc")]
  #[trace]
  fn test_truncate_with_custom_ellipsis(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] ellipsis: &str,
    #[case] direction: TruncateDirection,
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
  fn test_name_list(#[case] names: Vec<&str>, #[case] expected: &str) {
    let names = names.into_iter().map(str::to_owned).collect::<Vec<_>>();
    assert_eq!(expected, name_list(&names));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero( 0, "0 tips")]
  #[case::one(  1, "1 tip")]
  #[case::many( 2, "2 tips")]
  #[trace]
  fn test_counted(#[case] count: usize, #[case] expected: &str) {
    assert_eq!(expected, counted(count, "tip"));
  }
}
