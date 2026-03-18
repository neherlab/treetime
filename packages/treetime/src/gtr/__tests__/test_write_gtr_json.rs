#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{GtrModelName, JC69Params, jc69, write_gtr_json};
  use rstest::rstest;
  use tempfile::TempDir;

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

    assert!(dir.path().join("gtr_sparse.json").exists(), "gtr_sparse.json should exist");
    assert!(dir.path().join("gtr_dense.json").exists(), "gtr_dense.json should exist");
    assert!(
      !dir.path().join("gtr.json").exists(),
      "gtr.json should not exist when qualifiers are used"
    );
  }
}
