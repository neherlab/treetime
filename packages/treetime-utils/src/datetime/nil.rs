use regex::Regex;
use std::sync::LazyLock;

pub fn is_nil(input: &str) -> bool {
  const NON_VALUES: &[&str] = &[
    "na",
    "n/a",
    "nan",
    "null",
    "nil",
    "none",
    "empty",
    "missing",
    "undefined",
  ];
  #[allow(clippy::unwrap_used, reason = "compile-time-constant pattern; a malformed literal is a build-time bug")]
  static REGEX: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\?+|-+)$").unwrap());
  let input = input.to_lowercase();
  input.is_empty() || NON_VALUES.contains(&input.as_str()) || REGEX.is_match(&input)
}
