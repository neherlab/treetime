#[cfg(test)]
mod tests {
  use crate::io::json::{JsonPretty, json_read_file, json_write_file};
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
}
