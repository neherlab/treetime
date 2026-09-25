#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::runs::inputs::{hash_inputs, sha256_hex};
  use crate::runs::record::RunInput;
  use helpers::settings;
  use pretty_assertions::{assert_eq, assert_ne};
  use serde_json::{Map, json};
  use std::fs;
  use std::path::Path;
  use tempfile::tempdir;

  #[test]
  fn test_inputs_sha256_of_known_text() {
    assert_eq!(
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
      sha256_hex(b"abc").unwrap()
    );
  }

  #[test]
  fn test_inputs_record_path_size_and_sha256() {
    let dir = tempdir().unwrap();
    let tree = dir.path().join("tree.nwk");
    fs::write(&tree, "abc").unwrap();
    let hashed = hash_inputs(AppCommand::Clock, &settings(json!({ "tree": tree }))).unwrap();
    assert_eq!(
      vec![RunInput {
        setting: "tree".to_owned(),
        path: tree,
        size: 3,
        sha256: "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad".to_owned(),
      }],
      hashed.inputs
    );
  }

  #[test]
  fn test_inputs_config_hash_ignores_input_locations_and_outputs() {
    let first = tempdir().unwrap();
    let second = tempdir().unwrap();
    for dir in [first.path(), second.path()] {
      fs::write(dir.join("tree.nwk"), "(A:1,B:1);\n").unwrap();
      fs::write(dir.join("metadata.tsv"), "name\tdate\nA\t2000\nB\t2001\n").unwrap();
    }
    let config = |dir: &Path, out: &str| {
      settings(json!({
        "tree": dir.join("tree.nwk"),
        "metadata": dir.join("metadata.tsv"),
        "clock_filter": 3.0,
        "output_all": out,
      }))
    };
    let a = hash_inputs(AppCommand::Clock, &config(first.path(), "/runs/a/out")).unwrap();
    let b = hash_inputs(AppCommand::Clock, &config(second.path(), "/runs/b/out")).unwrap();
    assert_eq!(a.config_hash, b.config_hash);
  }

  #[test]
  fn test_inputs_config_hash_changes_with_input_contents_and_settings() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("a.nwk"), "(A:1,B:1);\n").unwrap();
    fs::write(dir.path().join("b.nwk"), "(A:2,B:1);\n").unwrap();
    let hash = |tree: &str, clock_filter: f64| {
      hash_inputs(
        AppCommand::Clock,
        &settings(json!({ "tree": dir.path().join(tree), "clock_filter": clock_filter })),
      )
      .unwrap()
      .config_hash
    };
    let base = hash("a.nwk", 3.0);
    assert_ne!(base, hash("b.nwk", 3.0));
    assert_ne!(base, hash("a.nwk", 2.0));
    assert_eq!(base, hash("a.nwk", 3.0));
  }

  #[test]
  fn test_inputs_config_hash_does_not_depend_on_key_order() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("t.nwk"), "(A:1,B:1);\n").unwrap();
    let tree = dir.path().join("t.nwk");
    let forward = settings(json!({ "tree": tree, "clock_filter": 3.0, "keep_root": true }));
    let mut backward = Map::new();
    for (key, value) in forward.iter().rev() {
      backward.insert(key.clone(), value.clone());
    }
    assert_eq!(
      hash_inputs(AppCommand::Clock, &forward).unwrap().config_hash,
      hash_inputs(AppCommand::Clock, &backward).unwrap().config_hash
    );
  }

  mod helpers {
    use serde_json::{Map, Value};

    pub(super) fn settings(value: Value) -> Map<String, Value> {
      match value {
        Value::Object(settings) => settings,
        other => panic!("expected a mapping of settings, got {other}"),
      }
    }
  }
}
