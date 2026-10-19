use itertools::Itertools;
use std::fmt::Display;

const NAMES_SHOWN: usize = 3;

#[macro_export]
macro_rules! o {
  ($x:expr $(,)?) => {
    ToOwned::to_owned($x)
  };
}

pub fn quote(x: impl Display) -> String {
  format!("\"{x}\"")
}

pub fn quote_single(x: impl Display) -> String {
  format!("'{x}'")
}

pub fn name_list(names: &[String]) -> String {
  let shown = names.iter().take(NAMES_SHOWN).join(", ");
  if names.len() > NAMES_SHOWN {
    format!("{shown}, and {} more", names.len() - NAMES_SHOWN)
  } else {
    shown
  }
}

pub fn counted(count: usize, noun: &str) -> String {
  if count == 1 {
    format!("{count} {noun}")
  } else {
    format!("{count} {noun}s")
  }
}

pub fn truncate_right_with_ellipsis(s: impl AsRef<str>, max_len: usize) -> String {
  truncate(s, max_len, Some("..."), TruncateDirection::Right)
}

#[expect(
  clippy::string_slice,
  reason = "the input is ASCII, so every byte index is a char boundary"
)]
#[expect(
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TruncateDirection {
  Left,
  Right,
  Middle,
}
