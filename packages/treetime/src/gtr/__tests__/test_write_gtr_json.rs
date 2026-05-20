#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{GtrModelName, GtrOutput, JC69Params, jc69, write_gtr_json};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::TempDir;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  #[rustfmt::skip]
  #[rstest]
  #[case::no_qualifier(   None,             "gtr.json")]
  #[case::sparse(         Some("sparse"),   "gtr_sparse.json")]
  #[case::dense(          Some("dense"),    "gtr_dense.json")]
  #[trace]
  fn test_write_gtr_json_filename(#[case] qualifier: Option<&str>, #[case] expected_filename: &str) {
    let dir = TempDir::new().unwrap();
    let gtr = jc69(JC69Params::default()).unwrap();
    write_gtr_json(&gtr, GtrModelName::JC69, dir.path(), qualifier).unwrap();

    let expected_path = dir.path().join(expected_filename);
    assert!(expected_path.exists(), "Expected file {expected_filename} not found");
  }

  #[test]
  fn test_write_gtr_json_both_partitions_no_overwrite() {
    let dir = TempDir::new().unwrap();
    let gtr = jc69(JC69Params::default()).unwrap();

    write_gtr_json(&gtr, GtrModelName::JC69, dir.path(), Some("sparse")).unwrap();
    write_gtr_json(&gtr, GtrModelName::JC69, dir.path(), Some("dense")).unwrap();

    assert!(
      dir.path().join("gtr_sparse.json").exists(),
      "gtr_sparse.json should exist"
    );
    assert!(
      dir.path().join("gtr_dense.json").exists(),
      "gtr_dense.json should exist"
    );
    assert!(
      !dir.path().join("gtr.json").exists(),
      "gtr.json should not exist when qualifiers are used"
    );
  }

  #[test]
  fn test_gtr_output_without_discrete_states_omits_fields() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let output = GtrOutput::new(&gtr, GtrModelName::JC69);
    let json = json_write_str(&output, JsonPretty(false))?;
    let parsed: serde_json::Value = serde_json::from_str(&json)?;

    assert!(parsed.get("attribute").is_none(), "attribute should be absent when None");
    assert!(parsed.get("states").is_none(), "states should be absent when None");
    assert_eq!(parsed["n_states"], 4);
    assert_eq!(parsed["model_name"], "JC69");
    assert_eq!(parsed["model_type"], "named");

    Ok(())
  }

  #[test]
  fn test_gtr_output_with_discrete_states_includes_fields() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let output = GtrOutput::new(&gtr, GtrModelName::Infer)
      .with_discrete_states("country", ["france", "germany", "usa"].iter());
    let json = json_write_str(&output, JsonPretty(false))?;
    let parsed: serde_json::Value = serde_json::from_str(&json)?;

    assert_eq!(parsed["attribute"], "country");
    assert_eq!(
      parsed["states"],
      serde_json::json!(["france", "germany", "usa"])
    );
    assert_eq!(parsed["n_states"], 4);
    assert_eq!(parsed["model_type"], "inferred");

    Ok(())
  }

  #[test]
  fn test_gtr_output_with_discrete_states_roundtrip() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let original = GtrOutput::new(&gtr, GtrModelName::Infer)
      .with_discrete_states("region", ["asia", "europe"].iter());
    let json = json_write_str(&original, JsonPretty(true))?;
    let restored: GtrOutput = serde_json::from_str(&json)?;

    assert_eq!(restored.attribute, Some("region".to_owned()));
    assert_eq!(
      restored.states,
      Some(vec!["asia".to_owned(), "europe".to_owned()])
    );
    assert_eq!(restored.model_type, original.model_type);
    assert_eq!(restored.model_name, original.model_name);
    assert_eq!(restored.n_states, original.n_states);

    Ok(())
  }
}
