#[cfg(test)]
mod tests {
  use crate::io::json::{
    SerializationFormat, json_or_yaml_read_file, json_or_yaml_write_file, serialization_format_from_path,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde::{Deserialize, Serialize};
  use tempfile::tempdir;

  #[rustfmt::skip]
  #[rstest]
  #[case::json(           "result.json",    SerializationFormat::Json)]
  #[case::json_compressed("result.json.xz", SerializationFormat::Json)]
  #[case::yaml(           "result.yaml",    SerializationFormat::Yaml)]
  #[case::yaml_compressed("result.yaml.xz", SerializationFormat::Yaml)]
  #[case::yml_compressed( "result.yml.gz",  SerializationFormat::Yaml)]
  #[trace]
  fn test_json_serialization_format_from_path(
    #[case] filepath: &str,
    #[case] expected: SerializationFormat,
  ) {
    let actual = serialization_format_from_path(filepath);
    assert_eq!(expected, actual);
  }

  #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
  struct Sample {
    name: String,
    count: u32,
  }

  // `json_or_yaml_read_file` must round-trip whatever `json_or_yaml_write_file` produced, across both
  // formats and through the transparent (de)compression that the underlying file helpers apply.
  #[rustfmt::skip]
  #[rstest]
  #[case::json(          "config.json")]
  #[case::json_xz(       "config.json.xz")]
  #[case::yaml(          "config.yaml")]
  #[case::yml_gz(        "config.yml.gz")]
  #[trace]
  fn test_json_or_yaml_read_file_roundtrip(#[case] filename: &str) {
    let dir = tempdir().unwrap();
    let path = dir.path().join(filename);
    let expected = Sample { name: "flu".to_owned(), count: 42 };
    json_or_yaml_write_file(&path, &expected).unwrap();
    let actual: Sample = json_or_yaml_read_file(&path).unwrap();
    assert_eq!(expected, actual);
  }
}
