#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{GtrModelName, JC69Params, jc69, write_gtr_json};
  use rstest::rstest;
  use std::sync::atomic::{AtomicU64, Ordering};

  static COUNTER: AtomicU64 = AtomicU64::new(0);

  fn unique_test_dir() -> std::path::PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("treetime_test_write_gtr_{id}"));
    std::fs::create_dir_all(&dir).unwrap();
    dir
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_qualifier(   None,             "gtr.json")]
  #[case::sparse(         Some("sparse"),   "gtr_sparse.json")]
  #[case::dense(          Some("dense"),    "gtr_dense.json")]
  #[trace]
  fn test_write_gtr_json_filename(#[case] qualifier: Option<&str>, #[case] expected_filename: &str) {
    let dir = unique_test_dir();
    let gtr = jc69(JC69Params::default()).unwrap();
    write_gtr_json(&gtr, GtrModelName::JC69, &dir, qualifier).unwrap();

    let expected_path = dir.join(expected_filename);
    assert!(expected_path.exists(), "Expected file {expected_filename} not found");

    std::fs::remove_dir_all(&dir).unwrap();
  }

  #[test]
  fn test_write_gtr_json_both_partitions_no_overwrite() {
    let dir = unique_test_dir();
    let gtr = jc69(JC69Params::default()).unwrap();

    write_gtr_json(&gtr, GtrModelName::JC69, &dir, Some("sparse")).unwrap();
    write_gtr_json(&gtr, GtrModelName::JC69, &dir, Some("dense")).unwrap();

    assert!(dir.join("gtr_sparse.json").exists(), "gtr_sparse.json should exist");
    assert!(dir.join("gtr_dense.json").exists(), "gtr_dense.json should exist");
    assert!(
      !dir.join("gtr.json").exists(),
      "gtr.json should not exist when qualifiers are used"
    );

    std::fs::remove_dir_all(&dir).unwrap();
  }
}
