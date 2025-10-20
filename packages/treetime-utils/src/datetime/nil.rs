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
  static REGEX: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\?+|-+)$").unwrap());
  let input = input.to_lowercase();
  input.is_empty() || NON_VALUES.contains(&input.as_str()) || REGEX.is_match(&input)
}
