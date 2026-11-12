#[cfg(test)]
mod tests {
  use app_backend::runs::manager::RunManager;
  use app_backend::runs::record::{CreateRunRequest, RunStatus, StartRunRequest};
  use app_commands::command::AppCommand;
  use app_commands::json_value::SparseConfig;
  use helpers::{extract, output_texts, run_cli_line, upload};
  use pretty_assertions::assert_eq;
  use serde_json::{Map, Value, json};
  use std::fs;
  use tempfile::tempdir;

  #[test]
  fn test_run_package_cli_runs_the_yaml_of_an_extracted_package_and_writes_equal_outputs() {
    let work = tempdir().unwrap();
    let runs = RunManager::open(&work.path().join("runs")).unwrap();
    let created = runs
      .create(&CreateRunRequest {
        command: AppCommand::Timetree,
        config: SparseConfig(Map::new()),
        defer_start: true,
        source_run_id: None,
      })
      .unwrap();
    let mut config = Map::new();
    config.insert("max_iter".to_owned(), json!(1));
    config.insert("tree".to_owned(), upload(&runs, &created.id, "tree.nwk"));
    config.insert("metadata".to_owned(), upload(&runs, &created.id, "metadata.tsv"));
    config.insert(
      "alignment".to_owned(),
      Value::Array(vec![upload(&runs, &created.id, "aln.fasta.xz")]),
    );
    let start = StartRunRequest {
      command: None,
      config: Some(SparseConfig(config)),
      source_run_id: None,
    };
    runs.start(&created.id, &start, Box::new(|_| Ok(()))).unwrap().run();
    let record = runs.get(&created.id).unwrap();
    let mut bytes = vec![];
    runs.package(&record.id).unwrap().stream(&mut bytes).unwrap();
    let folder = extract(&bytes, &work.path().join("package")).join(format!("treetime-{}", record.id.as_str()));
    let stored = record
      .inputs
      .iter()
      .filter_map(|input| {
        let key = input.stored()?;
        Some((
          runs.input_store().path(&key).to_string_lossy().into_owned(),
          input.path.to_string_lossy().into_owned(),
        ))
      })
      .collect::<Vec<_>>();
    let packaged = output_texts(&folder.join("out"), &stored);
    fs::remove_dir_all(folder.join("out")).unwrap();

    run_cli_line(
      &folder,
      fs::read_to_string(folder.join("command.sh"))
        .unwrap()
        .lines()
        .last()
        .unwrap(),
    );

    let extracted = vec![(format!("{}/", folder.to_string_lossy()), String::new())];
    assert_eq!(
      (RunStatus::Ok, true, packaged),
      (
        record.status,
        record.seed.is_some(),
        output_texts(&folder.join("out"), &extracted)
      )
    );
  }

  mod helpers {
    use crate::cli::treetime_cli::{resolve_args, treetime_parse_cli_args};
    use crate::run::run_command;
    use app_backend::job::JobId;
    use app_backend::runs::manager::RunManager;
    use clap::ColorChoice;
    use serde_json::Value;
    use std::collections::BTreeMap;
    use std::fs;
    use std::io::Cursor;
    use std::path::{Path, PathBuf};
    use treetime::progress::NoopProgress;
    use zip::ZipArchive;

    pub(super) fn upload(runs: &RunManager, id: &JobId, name: &str) -> Value {
      let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      let bytes = fs::read(data.join(name)).unwrap();
      let uploaded = runs.upload_input(id, name, &[], &mut &bytes[..], usize::MAX).unwrap();
      Value::String(uploaded.path.to_string_lossy().into_owned())
    }

    pub(super) fn extract(bytes: &[u8], target: &Path) -> PathBuf {
      ZipArchive::new(Cursor::new(bytes)).unwrap().extract(target).unwrap();
      target.to_path_buf()
    }

    pub(super) fn run_cli_line(folder: &Path, line: &str) {
      let argv = line
        .split_whitespace()
        .map(|token| {
          let path = folder.join(token);
          if path.is_file() {
            path.to_string_lossy().into_owned()
          } else {
            token.to_owned()
          }
        })
        .collect::<Vec<_>>();
      let parsed = treetime_parse_cli_args(argv, ColorChoice::Never).unwrap();
      let args = resolve_args(parsed).unwrap();
      run_command(args.command, &NoopProgress, &NoopProgress).unwrap();
    }

    pub(super) fn output_texts(dir: &Path, replacements: &[(String, String)]) -> BTreeMap<String, String> {
      fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_file())
        .map(|path| {
          let name = path.file_name().unwrap().to_string_lossy().into_owned();
          let text = replacements
            .iter()
            .fold(normalized_content(&path), |text, (from, to)| text.replace(from, to));
          (name, text)
        })
        .collect()
    }

    fn normalized_content(path: &Path) -> String {
      let name = path.file_name().unwrap().to_string_lossy();
      let bytes = fs::read(path).unwrap();
      if name.ends_with(".auspice.json") {
        let mut value: Value = serde_json::from_slice(&bytes).unwrap();
        if let Some(meta) = value.get_mut("meta").and_then(Value::as_object_mut) {
          meta.remove("updated");
        }
        return value.to_string();
      }
      String::from_utf8_lossy(&bytes).into_owned()
    }
  }
}
