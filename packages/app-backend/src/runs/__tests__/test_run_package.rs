#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{sparse, unconfined};
  use crate::results::run_results::run_results;
  use crate::runs::__tests__::test_run_storage::tests::helpers::{ZIKA, run_with_uploads, start_from};
  use crate::runs::app_events::AppChange;
  use crate::runs::files::RUN_DOWNLOAD_CHUNK_SIZE;
  use crate::runs::manager::RunManager;
  use crate::runs::record::{CreateRunRequest, RunStatus, StartRunRequest};
  use app_commands::command::AppCommand;
  use helpers::{
    Damage, damaged, data_dir, entry_text, package_bytes, package_entries, record_json, staged_folders,
    subscribe_changes,
  };
  use indoc::indoc;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use std::fs;
  use std::io::Cursor;
  use tempfile::tempdir;
  use treetime_utils::assert_error;
  use zip::CompressionMethod;

  #[test]
  fn test_run_package_import_in_another_runs_folder_recreates_the_run_and_its_results() {
    let (source_root, target_root) = (tempdir().unwrap(), tempdir().unwrap());
    let source = RunManager::open(source_root.path()).unwrap();
    let record = run_with_uploads(
      &source,
      AppCommand::Timetree,
      ZIKA,
      &json!({ "max_iter": 1, "sequence_length": 10807 }),
    );
    let target = RunManager::open(target_root.path()).unwrap();

    let imported = target.import(&mut &package_bytes(&source, &record.id)[..]).unwrap();

    assert_eq!(
      (
        RunStatus::Ok,
        record_json(&source.get(&record.id).unwrap()),
        run_results(&source, &record.id).unwrap()
      ),
      (
        imported.status,
        record_json(&imported),
        run_results(&target, &imported.id).unwrap()
      )
    );
  }

  #[test]
  fn test_run_package_import_beside_its_source_gets_a_new_id_and_emits_one_run_created_event() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let bytes = package_bytes(&runs, &record.id);
    let changes = subscribe_changes(&runs);

    let imported = runs.import(&mut &bytes[..]).unwrap();

    let mut expected = record_json(&runs.get(&record.id).unwrap());
    expected["id"] = json!(imported.id.as_str());
    expected["group_id"] = json!(imported.id.as_str());
    let changes = changes.lock();
    let [AppChange::RunCreated { run }] = changes.as_slice() else {
      panic!("expected one run-created change, got {changes:?}");
    };
    assert_eq!(
      (
        expected,
        run_results(&runs, &record.id).unwrap(),
        imported.id.clone(),
        true
      ),
      (
        record_json(&imported),
        run_results(&runs, &imported.id).unwrap(),
        run.id.clone(),
        imported.id != record.id
      )
    );
  }

  #[test]
  fn test_run_package_holds_the_run_folder_the_inputs_the_yaml_config_and_its_command() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let folder = format!("treetime-{}", record.id.as_str());

    let entries = package_entries(&package_bytes(&runs, &record.id));

    let expected_names = ["run.json", "run.yaml", "command.sh", "events.jsonl"]
      .into_iter()
      .map(str::to_owned)
      .chain(["inputs/metadata.tsv".to_owned(), "inputs/tree.nwk".to_owned()])
      .chain(
        runs
          .files(&record.id)
          .unwrap()
          .into_iter()
          .map(|file| format!("out/{}", file.path)),
      )
      .map(|name| format!("{folder}/{name}"))
      .sorted()
      .collect_vec();
    let yaml = entry_text(&entries, &format!("{folder}/run.yaml"));
    assert_eq!(
      (
        expected_names,
        indoc! {r#"
          #!/bin/sh
          # Reproduces the run with the TreeTime command-line tool. Run it in this folder: the outputs replace the files in out/.
          treetime clock --config run.yaml
        "#}
        .to_owned(),
        (true, true, true),
      ),
      (
        entries.iter().map(|(name, _, _)| name.clone()).sorted().collect_vec(),
        entry_text(&entries, &format!("{folder}/command.sh")),
        (
          yaml.lines().any(|line| line == r#"tree: "inputs/tree.nwk""#),
          yaml.lines().any(|line| line == r#"metadata: "inputs/metadata.tsv""#),
          yaml.lines().any(|line| line == r#"output_all: "out""#),
        ),
      )
    );
  }

  #[test]
  fn test_run_package_copies_inputs_from_outside_the_runs_folder_and_the_imported_run_reruns() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(&root.path().join("desktop")).unwrap();
    let data = data_dir().join(ZIKA);
    let created = runs
      .create(&CreateRunRequest {
        command: AppCommand::Clock,
        config: sparse(json!({ "tree": data.join("tree.nwk"), "metadata": data.join("metadata.tsv") })),
        defer_start: true,
        source_run_id: None,
        example: None,
      })
      .unwrap();
    runs
      .start(&created.id, &StartRunRequest::default(), unconfined())
      .unwrap()
      .run();
    let record = runs.get(&created.id).unwrap();
    let bytes = package_bytes(&runs, &record.id);
    let packaged = entry_text(
      &package_entries(&bytes),
      &format!("treetime-{}/run.json", record.id.as_str()),
    );
    let packaged: serde_json::Value = serde_json::from_str(&packaged).unwrap();
    let target = RunManager::open(&root.path().join("web")).unwrap();

    let imported = target.import(&mut &bytes[..]).unwrap();
    let rerun = start_from(&target, &target.settings(&imported.id).unwrap(), Some(&imported.id));

    assert_eq!(
      (
        (RunStatus::Ok, true),
        (
          json!("inputs/tree.nwk"),
          json!("inputs/metadata.tsv"),
          vec![json!("inputs/metadata.tsv"), json!("inputs/tree.nwk")]
        ),
        (RunStatus::Ok, true),
      ),
      (
        (
          record.status,
          record.inputs.iter().all(|input| input.path.is_absolute())
        ),
        (
          packaged["config"]["tree"].clone(),
          packaged["config"]["metadata"].clone(),
          packaged["inputs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|input| input["path"].clone())
            .sorted_by_key(ToString::to_string)
            .collect_vec()
        ),
        (
          rerun.status,
          run_results(&target, &imported.id).unwrap().problems.is_empty()
        ),
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::entry_outside_the_package(Damage::EntryOutside)]
  #[case::input_path_outside(       Damage::InputPathOutside)]
  #[case::absolute_setting_path(    Damage::AbsoluteSetting)]
  #[case::input_size(               Damage::InputSize)]
  #[case::input_sha256(             Damage::InputSha256)]
  #[case::missing_output_file(      Damage::MissingOutput)]
  #[case::run_not_ended(            Damage::NotEnded)]
  #[case::not_a_zip_archive(        Damage::NotZip)]
  #[trace]
  fn test_run_package_import_refuses_a_damaged_package_and_leaves_no_partial_run(#[case] damage: Damage) {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let inputs = runs.input_store().keys().unwrap();
    let (bytes, message) = damaged(&package_bytes(&runs, &record.id), &record, damage);

    assert_error!(runs.import(&mut &bytes[..]), message);
    assert_eq!(
      (vec![record.id], inputs, 0),
      (
        runs.store().list().unwrap().into_iter().map(|run| run.id).collect_vec(),
        runs.input_store().keys().unwrap(),
        staged_folders(root.path()),
      )
    );
  }

  #[test]
  fn test_run_package_of_a_run_that_has_not_ended_is_refused() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = runs
      .create(&CreateRunRequest {
        command: AppCommand::Clock,
        config: sparse(json!({})),
        defer_start: true,
        source_run_id: None,
        example: None,
      })
      .unwrap();
    assert_error!(
      runs.package(&created.id).map(|_| ()),
      format!(
        "run `{}` has not ended (status `created`); a run package holds a run that has ended",
        created.id.as_str()
      )
    );
  }

  #[test]
  fn test_run_package_chunks_stay_within_the_chunk_size_and_store_every_entry() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let large = b"ACGT"
      .iter()
      .copied()
      .cycle()
      .take(3 * RUN_DOWNLOAD_CHUNK_SIZE + 17)
      .collect_vec();
    fs::write(runs.store().out_dir(&record.id).join("large.fasta"), &large).unwrap();

    let chunks: Vec<Vec<u8>> = runs.package(&record.id).unwrap().chunks().try_collect().unwrap();

    let entries = package_entries(&chunks.concat());
    let large_name = format!("treetime-{}/out/large.fasta", record.id.as_str());
    assert_eq!(
      (true, true, Some(large)),
      (
        chunks.iter().all(|chunk| chunk.len() <= RUN_DOWNLOAD_CHUNK_SIZE),
        entries
          .iter()
          .all(|(_, method, _)| *method == CompressionMethod::Stored),
        entries
          .iter()
          .find(|(name, _, _)| *name == large_name)
          .map(|(_, _, content)| content.clone()),
      )
    );
  }

  #[test]
  fn test_run_package_written_to_a_file_deflates_the_entries_of_the_stream() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let mut written = Cursor::new(vec![]);

    runs.package(&record.id).unwrap().write(&mut written).unwrap();

    let streamed = package_entries(&package_bytes(&runs, &record.id));
    let written = package_entries(&written.into_inner());
    assert_eq!(
      (
        streamed
          .iter()
          .map(|(name, _, content)| (name.clone(), content.clone()))
          .collect_vec(),
        true
      ),
      (
        written
          .iter()
          .map(|(name, _, content)| (name.clone(), content.clone()))
          .collect_vec(),
        written
          .iter()
          .all(|(_, method, _)| *method == CompressionMethod::Deflated)
      )
    );
  }

  #[test]
  fn test_run_package_chunks_end_after_a_file_that_cannot_be_read() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let mut chunks = runs.package(&record.id).unwrap().chunks();
    let removed = runs.store().out_dir(&record.id).join("clock.clock-model.json");
    fs::remove_file(&removed).unwrap();

    let mut items = chunks.by_ref().collect_vec();
    let last = items.pop().unwrap();

    assert_eq!(
      (true, None),
      (
        items.iter().all(Result::is_ok),
        chunks.next().map(|chunk| chunk.is_ok())
      )
    );
    assert_error!(
      last,
      format!(
        "When opening '{}': No such file or directory (os error 2)",
        removed.display()
      )
    );
  }

  mod helpers {
    use crate::job::JobId;
    use crate::runs::app_events::AppChange;
    use crate::runs::manager::RunManager;
    use crate::runs::record::RunRecord;
    use parking_lot::Mutex;
    use serde_json::{Value, json};
    use std::fs;
    use std::io::{Cursor, Read, Write};
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use treetime_utils::io::json::to_json_value;
    use zip::write::SimpleFileOptions;
    use zip::{CompressionMethod, ZipArchive, ZipWriter};

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Damage {
      EntryOutside,
      InputPathOutside,
      AbsoluteSetting,
      InputSize,
      InputSha256,
      MissingOutput,
      NotEnded,
      NotZip,
    }

    pub(super) fn damaged(bytes: &[u8], record: &RunRecord, damage: Damage) -> (Vec<u8>, String) {
      let folder = format!("treetime-{}", record.id.as_str());
      let tree = record
        .inputs
        .iter()
        .find(|input| input.path == Path::new("inputs/tree.nwk"))
        .unwrap();
      let outside = |path: &str| {
        format!(
          "the run package names `{path}`, a path outside the package; every path of a run package is relative and \
           inside it"
        )
      };
      let mismatch = |listed_size: usize, listed_sha256: &str| {
        format!(
          "input `inputs/tree.nwk` of the run package does not match its run.json: the file has {} bytes and SHA-256 {}, \
           run.json lists {listed_size} bytes and SHA-256 {listed_sha256}",
          tree.size, tree.sha256
        )
      };
      let edit_tree = |edit: &dyn Fn(&mut Value)| {
        edit_record(bytes, &folder, |record| {
          let input = record["inputs"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|input| input["path"] == json!("inputs/tree.nwk"))
            .unwrap();
          edit(input);
        })
      };
      match damage {
        Damage::EntryOutside => {
          let name = format!("{folder}/../evil.txt");
          (
            rewrite(bytes, |_, content| Some(content), &[(&name, b"evil")]),
            outside(&name),
          )
        },
        Damage::InputPathOutside => (
          edit_tree(&|input| input["path"] = json!("../tree.nwk")),
          outside("../tree.nwk"),
        ),
        Damage::AbsoluteSetting => (
          edit_record(bytes, &folder, |record| record["config"]["tree"] = json!("/etc/hosts")),
          outside("/etc/hosts"),
        ),
        Damage::InputSize => (
          edit_tree(&|input| input["size"] = json!(tree.size + 1)),
          mismatch(tree.size + 1, &tree.sha256),
        ),
        Damage::InputSha256 => {
          let sha256 = "0".repeat(64);
          (
            edit_tree(&|input| input["sha256"] = json!(sha256)),
            mismatch(tree.size, &sha256),
          )
        },
        Damage::MissingOutput => {
          let output = record.output_files[0].path.to_string_lossy().into_owned();
          let removed = format!("{folder}/out/{output}");
          (
            rewrite(bytes, |name, content| (name != removed).then_some(content), &[]),
            format!("the run package holds no output file `out/{output}`, which its run.json lists"),
          )
        },
        Damage::NotEnded => (
          edit_record(bytes, &folder, |record| record["status"] = json!("running")),
          format!(
            "the run package holds run `{}`, which has not ended (status `running`); a run package holds a run that \
             has ended",
            record.id.as_str()
          ),
        ),
        Damage::NotZip => (
          b"not a zip archive".to_vec(),
          "the file is not a run package: it is not a zip archive: invalid Zip archive: Could not find EOCD".to_owned(),
        ),
      }
    }

    fn edit_record(bytes: &[u8], folder: &str, edit: impl Fn(&mut Value)) -> Vec<u8> {
      let run_file = format!("{folder}/run.json");
      rewrite(
        bytes,
        |name, content| {
          if name != run_file {
            return Some(content);
          }
          let mut record: Value = serde_json::from_slice(&content).unwrap();
          edit(&mut record);
          Some(serde_json::to_vec_pretty(&record).unwrap())
        },
        &[],
      )
    }

    fn rewrite(bytes: &[u8], edit: impl Fn(&str, Vec<u8>) -> Option<Vec<u8>>, extra: &[(&str, &[u8])]) -> Vec<u8> {
      let mut archive = ZipArchive::new(Cursor::new(bytes)).unwrap();
      let mut writer = ZipWriter::new(Cursor::new(vec![]));
      for index in 0..archive.len() {
        let mut file = archive.by_index(index).unwrap();
        let name = file.name().unwrap().into_owned();
        let mut content = vec![];
        file.read_to_end(&mut content).unwrap();
        if let Some(content) = edit(&name, content) {
          writer.start_file(name, SimpleFileOptions::default()).unwrap();
          writer.write_all(&content).unwrap();
        }
      }
      for (name, content) in extra {
        writer.start_file(*name, SimpleFileOptions::default()).unwrap();
        writer.write_all(content).unwrap();
      }
      writer.finish().unwrap().into_inner()
    }

    pub(super) fn package_bytes(runs: &RunManager, id: &JobId) -> Vec<u8> {
      let mut bytes = vec![];
      runs.package(id).unwrap().stream(&mut bytes).unwrap();
      bytes
    }

    pub(super) fn package_entries(bytes: &[u8]) -> Vec<(String, CompressionMethod, Vec<u8>)> {
      let mut archive = ZipArchive::new(Cursor::new(bytes)).unwrap();
      (0..archive.len())
        .map(|index| {
          let mut file = archive.by_index(index).unwrap();
          let mut content = vec![];
          file.read_to_end(&mut content).unwrap();
          (file.name().unwrap().into_owned(), file.compression(), content)
        })
        .collect()
    }

    pub(super) fn entry_text(entries: &[(String, CompressionMethod, Vec<u8>)], name: &str) -> String {
      let (_, _, content) = entries.iter().find(|(entry, _, _)| entry == name).unwrap();
      String::from_utf8(content.clone()).unwrap()
    }

    pub(super) fn record_json(record: &RunRecord) -> Value {
      to_json_value(record).unwrap()
    }

    pub(super) fn subscribe_changes(runs: &RunManager) -> Arc<Mutex<Vec<AppChange>>> {
      let changes = Arc::new(Mutex::new(vec![]));
      let seen = Arc::clone(&changes);
      runs.app_events().subscribe(
        None,
        Box::new(move |event| {
          seen.lock().push(event.change.clone());
          true
        }),
      );
      changes
    }

    pub(super) fn staged_folders(root: &Path) -> usize {
      let staging = root.join(".staging");
      if staging.exists() {
        fs::read_dir(staging).unwrap().count()
      } else {
        0
      }
    }

    pub(super) fn data_dir() -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .canonicalize()
        .unwrap()
    }
  }
}
