use regex::regex;

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
  let input = input.to_lowercase();
  input.is_empty() || NON_VALUES.contains(&input.as_str()) || regex!(r"^(\?+|-+)$").is_match(&input)
}
