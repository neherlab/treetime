#[cfg(test)]
mod tests {
  use crate::runs::input_store::InputStore;
  use crate::runs::inputs::{InputHashes, hash_inputs, sha256_hex};
  use crate::runs::record::RunInput;
  use app_commands::command::AppCommand;
  use helpers::settings;
  #[cfg(unix)]
  use helpers::{Replacement, replace_keeping_modification_time, wait_for_a_later_change_time};
  use pretty_assertions::{assert_eq, assert_ne};
  use rstest::rstest;
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

  #[cfg(unix)]
  #[rustfmt::skip]
  #[rstest]
  #[case::rewritten_in_place(Replacement::InPlace)]
  #[case::renamed_into_place(Replacement::Renamed)]
  #[trace]
  fn test_inputs_digest_of_a_same_length_replacement_with_the_modification_time_restored(
    #[case] replacement: Replacement,
  ) {
    let dir = tempdir().unwrap();
    let store = InputStore::new(&dir.path().join("runs"));
    let tree = dir.path().join("tree.nwk");
    fs::write(&tree, "abc").unwrap();
    let hashes = InputHashes::default();
    let before = hashes.sha256(&store, &tree).unwrap();

    wait_for_a_later_change_time(dir.path(), &tree);
    replace_keeping_modification_time(&tree, b"xyz", replacement);

    assert_eq!(
      (sha256_hex(b"abc").unwrap(), sha256_hex(b"xyz").unwrap()),
      (before, hashes.sha256(&store, &tree).unwrap())
    );
  }

  mod helpers {
    use serde_json::{Map, Value};
    #[cfg(unix)]
    use std::fs::{self, OpenOptions};
    #[cfg(unix)]
    use std::io::Write;
    #[cfg(unix)]
    use std::os::unix::fs::MetadataExt;
    #[cfg(unix)]
    use std::path::Path;

    #[cfg(unix)]
    #[derive(Clone, Copy, Debug)]
    pub(super) enum Replacement {
      InPlace,
      Renamed,
    }

    #[cfg(unix)]
    pub(super) fn replace_keeping_modification_time(path: &Path, content: &[u8], replacement: Replacement) {
      let modified = fs::metadata(path).unwrap().modified().unwrap();
      let target = match replacement {
        Replacement::InPlace => path.to_path_buf(),
        Replacement::Renamed => path.with_extension("new"),
      };
      let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(&target)
        .unwrap();
      file.write_all(content).unwrap();
      file.set_modified(modified).unwrap();
      drop(file);
      if matches!(replacement, Replacement::Renamed) {
        fs::rename(&target, path).unwrap();
      }
      assert_eq!(modified, fs::metadata(path).unwrap().modified().unwrap());
    }

    #[cfg(unix)]
    pub(super) fn wait_for_a_later_change_time(dir: &Path, path: &Path) {
      let changed = change_time(path);
      let probe = dir.join("probe");
      fs::write(&probe, "x").unwrap();
      while change_time(&probe) <= changed {
        fs::write(&probe, "x").unwrap();
      }
    }

    #[cfg(unix)]
    fn change_time(path: &Path) -> (i64, i64) {
      let metadata = fs::metadata(path).unwrap();
      (metadata.ctime(), metadata.ctime_nsec())
    }

    pub(super) fn settings(value: Value) -> Map<String, Value> {
      match value {
        Value::Object(settings) => settings,
        other => panic!("expected a mapping of settings, got {other}"),
      }
    }
  }
}
