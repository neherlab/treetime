#[cfg(test)]
mod tests {
  use crate::io::json::{SerializationFormat, serialization_format_from_path};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

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
}
