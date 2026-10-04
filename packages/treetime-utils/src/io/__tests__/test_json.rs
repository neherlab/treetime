#[cfg(test)]
mod tests {
  use crate::assert_error;
  use crate::io::json::{JsonPretty, json_read_file, json_read_str, json_write_file, json_write_str};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde::{Deserialize, Serialize};
  use tempfile::tempdir;

  #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
  struct Sample {
    name: String,
    count: u32,
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::json(   "config.json")]
  #[case::json_xz("config.json.xz")]
  #[case::json_gz("config.json.gz")]
  #[trace]
  fn test_json_read_file_roundtrip(#[case] filename: &str) {
    let dir = tempdir().unwrap();
    let path = dir.path().join(filename);
    let expected = Sample { name: "flu".to_owned(), count: 42 };
    json_write_file(&path, &expected, JsonPretty(true)).unwrap();
    let actual: Sample = json_read_file(&path).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_json_read_str_rejects_trailing_characters() {
    assert_error!(
      json_read_str::<Sample>(r#"{"name":"flu","count":1} x"#),
      "When parsing JSON: trailing characters at line 1 column 26"
    );
  }

  #[test]
  fn test_json_write_str_returns_the_value_without_a_newline() -> Result<(), Report> {
    let sample = Sample {
      name: "flu".to_owned(),
      count: 1,
    };
    assert_eq!(
      r#"{"name":"flu","count":1}"#,
      json_write_str(&sample, JsonPretty(false))?
    );
    Ok(())
  }
}
