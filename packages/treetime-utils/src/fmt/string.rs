use std::fmt::Display;

#[macro_export]
macro_rules! o {
  ($x:expr $(,)?) => {
    ToOwned::to_owned($x)
  };
}

#[allow(
  clippy::as_conversions,
  clippy::unwrap_used,
  reason = "ASCII char sequence narrowed to bytes; the resulting buffer is valid UTF-8 by construction"
)]
pub fn vec_to_string(v: Vec<char>) -> String {
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

#[allow(clippy::string_slice)]
#[allow(
  clippy::integer_division,
  reason = "intentional halving of the budget for middle truncation; the remainder is assigned to the right half"
)]
pub fn truncate(s: impl AsRef<str>, max_len: usize, ellipsis: Option<&str>, direction: TruncateDirection) -> String {
  let s = s.as_ref();
  debug_assert!(s.is_ascii(), "Input to truncate must be ASCII");
  if s.len() <= max_len {
    return s.into();
  }

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

pub fn truncate_right(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, None, TruncateDirection::Right)
}

pub fn truncate_left(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, None, TruncateDirection::Left)
}

pub fn truncate_middle(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, None, TruncateDirection::Middle)
}

pub fn truncate_right_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Right)
}

pub fn truncate_left_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Left)
}

pub fn truncate_middle_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Middle)
}
