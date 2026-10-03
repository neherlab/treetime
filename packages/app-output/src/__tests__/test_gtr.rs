#[cfg(test)]
mod tests {
  use crate::gtr::write_gtr_json;
  use rstest::rstest;
  use tempfile::TempDir;
  use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, JC69Params, jc69};

  #[rustfmt::skip]
  #[rstest]
  #[case::no_qualifier(   "gtr.json")]
  #[case::sparse(         "gtr_sparse.json")]
  #[case::dense(          "gtr_dense.json")]
  #[trace]
  fn test_write_gtr_json_filename(#[case] filename: &str) {
    let dir = TempDir::new().unwrap();
    let gtr = jc69(JC69Params::default()).unwrap();
    let output = GtrOutput::builder().gtr(&gtr).model_name(GtrModelName::JC69).build();
    write_gtr_json(&output, dir.path().join(filename)).unwrap();

    let expected_path = dir.path().join(filename);
    assert!(expected_path.exists(), "Expected file {filename} not found");
  }

  #[test]
  fn test_write_gtr_json_both_partitions_no_overwrite() {
    let dir = TempDir::new().unwrap();
    let gtr = jc69(JC69Params::default()).unwrap();
    let output = GtrOutput::builder().gtr(&gtr).model_name(GtrModelName::JC69).build();

    write_gtr_json(&output, dir.path().join("gtr_sparse.json")).unwrap();
    write_gtr_json(&output, dir.path().join("gtr_dense.json")).unwrap();

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
}
