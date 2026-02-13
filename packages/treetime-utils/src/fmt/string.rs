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
