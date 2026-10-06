#[cfg(test)]
mod tests {
  use crate::config::source::{ConfigSource, parse_config_document};
  use eyre::Report;
  use helpers::{parse, parse_error_headline};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};

  #[test]
  fn test_source_parse_rejects_duplicate_mapping_key() {
    assert_eq!(
      "invalid configuration: could not parse config: Unexpected: duplicate map key \"a\" at line 2 column 1 (path: a)",
      parse_error_headline("a: 1\na: 2\n")
    );
  }

  #[test]
  fn test_source_parse_rejects_infinity() {
    assert_eq!(
      "invalid configuration: could not parse config: the number inf has no JSON representation at line 1 column 4",
      parse_error_headline("x: .inf\n")
    );
  }

  #[test]
  fn test_source_parse_rejects_nan() {
    assert_eq!(
      "invalid configuration: could not parse config: the number NaN has no JSON representation at line 1 column 4",
      parse_error_headline("x: .nan\n")
    );
  }

  #[test]
  fn test_source_parse_yaml11_booleans_no_and_on() {
    let value = parse("first: no\nsecond: on\n").unwrap();
    assert_eq!(json!({ "first": false, "second": true }), value);
  }

  #[test]
  fn test_source_parse_keeps_quoted_boolean_words_as_text() {
    let value = parse("first: \"no\"\nsecond: 'on'\n").unwrap();
    assert_eq!(json!({ "first": "no", "second": "on" }), value);
  }

  #[test]
  fn test_source_parse_applies_merge_key() {
    let value = parse("base: &anchor\n  shared: 1\nchild:\n  <<: *anchor\n  own: 2\n").unwrap();
    assert_eq!(json!({ "shared": 1 }), value["base"]);
    assert_eq!(json!({ "shared": 1, "own": 2 }), value["child"]);
  }

  #[test]
  fn test_source_parse_preserves_scientific_notation() {
    let value = parse("rate: 5.7e-05\n").unwrap();
    assert_eq!(json!({ "rate": 5.7e-05 }), value);
  }

  #[test]
  fn test_source_parse_reads_exponent_without_fraction_as_number() {
    let value = parse("prune_short: 1e-12\n").unwrap();
    assert_eq!(json!({ "prune_short": 1e-12 }), value);
  }

  mod helpers {
    use super::*;

    pub(super) fn parse(text: &str) -> Result<Value, Report> {
      let source = ConfigSource::new("config.yaml", text.to_owned());
      parse_config_document(&source, text)
    }

    pub(super) fn parse_error_headline(text: &str) -> String {
      let err = parse(text).expect_err("expected a parse error");
      err.to_string().lines().next().unwrap_or_default().to_owned()
    }
  }
}
