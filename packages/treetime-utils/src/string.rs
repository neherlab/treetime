use std::fmt::Display;

#[macro_export]
macro_rules! o {
  ($x:expr $(,)?) => {
    ToOwned::to_owned($x)
  };
}

pub fn vec_to_string(v: Vec<char>) -> String {
  // Surprisingly, this is the fastest way, according to `benches/vec_char_to_string.rs`
  let bytes: Vec<u8> = v.into_iter().map(|c| c as u8).collect();
  String::from_utf8(bytes).unwrap()
}

pub fn quote(x: impl Display) -> String {
  format!("\"{x}\"")
}

pub fn quote_single(x: impl Display) -> String {
  format!("'{x}'")
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TruncateDirection {
  Left,
  Right,
  Middle,
}

/// Truncates string to max_len bytes with ellipsis placement based on direction.
/// Assumes ASCII input for performance.
#[allow(clippy::string_slice)]
pub fn truncate(s: impl AsRef<str>, max_len: usize, ellipsis: Option<&str>, direction: TruncateDirection) -> String {
  let s = s.as_ref();
  debug_assert!(s.is_ascii(), "Input to truncate must be ASCII");
  if s.len() <= max_len {
    return s.into();
  }

  // If ellipsis is provided but too long, ignore it
  let ellipsis = ellipsis.filter(|ell| max_len >= ell.len());

  match (direction, ellipsis) {
    (TruncateDirection::Right, None) => s[..max_len].into(),
    (TruncateDirection::Right, Some(ell)) => [&s[..max_len - ell.len()], ell].concat(),
    (TruncateDirection::Left, None) => s[s.len() - max_len..].into(),
    (TruncateDirection::Left, Some(ell)) => [ell, &s[s.len() - (max_len - ell.len())..]].concat(),
    (TruncateDirection::Middle, None) => {
      let half = max_len / 2;
      [&s[..half], &s[s.len() - (max_len - half)..]].concat()
    },
    (TruncateDirection::Middle, Some(ell)) => {
      let remaining = max_len - ell.len();
      let left = remaining / 2;
      let right = remaining - left;
      [&s[..left], ell, &s[s.len() - right..]].concat()
    },
  }
}

/// Truncates string to max_len bytes from the right.
/// Assumes ASCII input for performance.
pub fn truncate_right(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, None, TruncateDirection::Right)
}

/// Truncates string to max_len bytes from the left.
/// Assumes ASCII input for performance.
pub fn truncate_left(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, None, TruncateDirection::Left)
}

/// Truncates string to max_len bytes from the middle.
/// Assumes ASCII input for performance.
pub fn truncate_middle(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, None, TruncateDirection::Middle)
}

/// Truncates string to max_len bytes from the right, appending "..." if truncated.
/// Assumes ASCII input for performance.
pub fn truncate_right_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Right)
}

/// Truncates string to max_len bytes from the left, prepending "..." if truncated.
/// Assumes ASCII input for performance.
pub fn truncate_left_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Left)
}

/// Truncates string to max_len bytes from the middle, inserting "..." if truncated.
/// Assumes ASCII input for performance.
pub fn truncate_middle_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Middle)
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case("hello world",                     8, "hello wo")]
  #[case("hello world",                    10, "hello worl")]
  #[case("hello world",                     4, "hell")]
  #[case("hello",                          10, "hello")]
  #[case("hello world",                     5, "hello")]
  #[case("hello world",                     1, "h")]
  #[case("",                                5, "")]
  #[case("abcdefghij",                      7, "abcdefg")]
  #[case("hello",                           0, "")]
  #[case("hello",                           5, "hello")]
  #[case("hello",                           6, "hello")]
  #[case("foo,bar;baz:qux",                10, "foo,bar;ba")]
  #[case("test@example.com",                8, "test@exa")]
  #[case("path/to/file.txt",               12, "path/to/file")]
  #[case("...already...has...ellipsis...", 15, "...already...ha")]
  #[case("!@#$%^&*()",                      5, "!@#$%")]
  #[case("a b c d e f g",                   7, "a b c d")]
  #[trace]
  fn test_truncate_right(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] expected: &str,
  ) {
    let actual = truncate_right(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("hello world",                     8, "lo world")]
  #[case("hello world",                    10, "ello world")]
  #[case("hello world",                     2, "ld")]
  #[case("hello",                          10, "hello")]
  #[case("hello world",                     5, "world")]
  #[case("hello world",                     1, "d")]
  #[case("",                                5, "")]
  #[case("abcdefghij",                      7, "defghij")]
  #[case("hello",                           0, "")]
  #[case("hello",                           5, "hello")]
  #[case("hello",                           6, "hello")]
  #[case("foo,bar;baz:qux",                10, "ar;baz:qux")]
  #[case("test@example.com",                8, "mple.com")]
  #[case("path/to/file.txt",               12, "/to/file.txt")]
  #[case("...already...has...ellipsis...", 15, "s...ellipsis...")]
  #[case("!@#$%^&*()",                      5, "^&*()")]
  #[case("a b c d e f g",                   7, "d e f g")]
  #[trace]
  fn test_truncate_left(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] expected: &str,
  ) {
    let actual = truncate_left(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("hello world",                     8, "hellorld")]
  #[case("hello world",                     4, "held")]
  #[case("hello world",                     9, "hellworld")]
  #[case("hello",                          10, "hello")]
  #[case("hello world",                     5, "herld")]
  #[case("hello world",                     1, "d")]
  #[case("",                                5, "")]
  #[case("abcdefghij",                      7, "abcghij")]
  #[case("hello",                           0, "")]
  #[case("hello",                           2, "ho")]
  #[case("hello",                           5, "hello")]
  #[case("hello",                           6, "hello")]
  #[case("foo,bar;baz:qux",                10, "foo,bz:qux")]
  #[case("test@example.com",                8, "test.com")]
  #[case("path/to/file.txt",               12, "path/tle.txt")]
  #[case("...already...has...ellipsis...", 15, "...alreipsis...")]
  #[case("!@#$%^&*()",                      5, "!@*()")]
  #[case("a b c d e f g",                   7, "a b f g")]
  #[trace]
  fn test_truncate_middle(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] expected: &str,
  ) {
    let actual = truncate_middle(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("hello world",                     8, "hello...")]
  #[case("hello world",                     4, "h...")]
  #[case("hello world",                    11, "hello world")]
  #[case("hello",                          10, "hello")]
  #[case("hello world",                     5, "he...")]
  #[case("hello world",                     3, "...")]
  #[case("",                                5, "")]
  #[case("abcdefghij",                      7, "abcd...")]
  #[case("hello",                           0, "")]
  #[case("hello",                           2, "he")]
  #[case("hello",                           1, "h")]
  #[case("hello",                           5, "hello")]
  #[case("hello",                           6, "hello")]
  #[case("foo,bar;baz:qux",                10, "foo,bar...")]
  #[case("test@example.com",               12, "test@exam...")]
  #[case("path/to/file.txt",               12, "path/to/f...")]
  #[case("...already...has...ellipsis...", 20, "...already...has....")]
  #[case("!@#$%^&*()",                      7, "!@#$...")]
  #[case("a b c d e f g",                  10, "a b c d...")]
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
  #[case("hello world",                     8, "...world")]
  #[case("hello world",                     4, "...d")]
  #[case("hello world",                    11, "hello world")]
  #[case("hello",                          10, "hello")]
  #[case("hello world",                     5, "...ld")]
  #[case("hello world",                     3, "...")]
  #[case("",                                5, "")]
  #[case("abcdefghij",                      7, "...ghij")]
  #[case("hello",                           0, "")]
  #[case("hello",                           2, "lo")]
  #[case("hello",                           1, "o")]
  #[case("hello",                           5, "hello")]
  #[case("hello",                           6, "hello")]
  #[case("foo,bar;baz:qux",                10, "...baz:qux")]
  #[case("test@example.com",               12, "...ample.com")]
  #[case("path/to/file.txt",               12, ".../file.txt")]
  #[case("...already...has...ellipsis...", 20, "...has...ellipsis...")]
  #[case("...starts...with...dots...",     18, "....with...dots...")]
  #[case("ends...with...dots...",          15, "...th...dots...")]
  #[case("..........many...dots..........", 20, "......dots..........")]
  #[case("!@#$%^&*()",                      7, "...&*()")]
  #[case("a b c d e f g",                  10, "...d e f g")]
  #[trace]
  fn test_truncate_left_with_ellipsis(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] expected: &str,
  ) {
    let actual = truncate_left_with_ellipsis(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("hello world",                     8, "he...rld")]
  #[case("hello world",                     9, "hel...rld")]
  #[case("hello world",                     4, "...d")]
  #[case("hello",                          10, "hello")]
  #[case("hello world",                     5, "h...d")]
  #[case("hello world",                     3, "...")]
  #[case("",                                5, "")]
  #[case("abcdefghij",                      7, "ab...ij")]
  #[case("hello",                           0, "")]
  #[case("hello",                           2, "ho")]
  #[case("hello",                           1, "o")]
  #[case("hello",                           5, "hello")]
  #[case("hello",                           6, "hello")]
  #[case("foo,bar;baz:qux",                10, "foo...:qux")]
  #[case("test@example.com",               12, "test...e.com")]
  #[case("path/to/file.txt",               12, "path...e.txt")]
  #[case("...already...has...ellipsis...", 20, "...alrea...lipsis...")]
  #[case("...starts...with...dots...",     18, "...star....dots...")]
  #[case("ends...with...dots...",          15, "ends.....ots...")]
  #[case("..........many...dots..........", 20, "....................")]
  #[case("!@#$%^&*()",                      7, "!@...()")]
  #[case("a b c d e f g",                  10, "a b... f g")]
  #[trace]
  fn test_truncate_middle_with_ellipsis(
    #[case] input: &str,
    #[case] max_len: usize,
    #[case] expected: &str,
  ) {
    let actual = truncate_middle_with_ellipsis(input, max_len);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("hello world",          8, "",       TruncateDirection::Right,  "hello wo")]
  #[case("hello world",          8, "<cut>",  TruncateDirection::Right,  "hel<cut>")]
  #[case("hello world",         10, "[...]",  TruncateDirection::Middle, "he[...]rld")]
  #[case("foo,bar;baz",          8, " >>",    TruncateDirection::Left,   " >>r;baz")]
  #[case("test@example.com",    12, ">>",     TruncateDirection::Left,   ">>xample.com")]
  #[case("path/to/file.txt",    15, "<snip>", TruncateDirection::Middle, "path<snip>e.txt")]
  #[case("a b c d e f g",       10, " ~ ",    TruncateDirection::Middle, "a b ~  f g")]
  #[case("!@#$%^&*()",           7, "***",    TruncateDirection::Right,  "!@#$***")]
  #[case("too short",           20, "...",    TruncateDirection::Right,  "too short")]
  #[case("exact",                5, "...",    TruncateDirection::Right,  "exact")]
  #[case("just over",            8, "....",   TruncateDirection::Middle, "ju....er")]
  #[case("ellipsis in input...", 15, "...",   TruncateDirection::Right,  "ellipsis in ...")]
  #[case("multiple words here", 12, " <> ",   TruncateDirection::Middle, "mult <> here")]
  #[case("x",                    1, "...",    TruncateDirection::Right,  "x")]
  #[case("xy",                   1, "...",    TruncateDirection::Left,   "y")]
  #[case("abc",                  2, "z",      TruncateDirection::Middle, "zc")]
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
}
