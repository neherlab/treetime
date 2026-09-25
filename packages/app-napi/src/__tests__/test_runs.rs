#[cfg(test)]
mod tests {
  use crate::runs::{create_run, start_run};
  use app_commands::job::TerminalEvent;
  use app_commands::runs::manager::RunManager;
  use app_commands::runs::record::RunStatus;
  use app_output::output_plan::OutputSelection;
  use helpers::ancestral_request;
  use pretty_assertions::assert_eq;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_runs_desktop_run_writes_its_outputs_into_the_run_folder() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = create_run(&runs, &ancestral_request()).unwrap();
    let started = start_run(&runs, record.id.as_str(), None).unwrap();
    let terminal = started.run();
    assert!(matches!(terminal, TerminalEvent::Ok { .. }), "{terminal:?}");

    let record = runs.get(&record.id).unwrap();
    let kinds: Vec<Option<OutputSelection>> = runs
      .files(&record.id)
      .unwrap()
      .into_iter()
      .filter(|file| file.path == "ancestral.auspice.json")
      .map(|file| file.kind)
      .collect();
    assert_eq!(
      (RunStatus::Ok, vec![Some(OutputSelection::Auspice)]),
      (record.status, kinds)
    );
  }

  #[test]
  fn test_runs_desktop_run_cannot_start_twice() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = create_run(&runs, &ancestral_request()).unwrap();
    let _started = start_run(&runs, record.id.as_str(), None).unwrap();
    assert_error!(
      start_run(&runs, record.id.as_str(), None),
      format!("run `{}` has already started", record.id.as_str())
    );
  }

  #[test]
  fn test_runs_desktop_rejects_malformed_ids() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    assert_error!(
      start_run(&runs, "../x", None),
      "invalid job id `../x`: expected 1 to 128 ASCII letters, digits, `-` or `_`"
    );
  }

  mod helpers {
    use serde_json::json;
    use std::path::Path;

    pub(super) fn ancestral_request() -> String {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({
        "command": "ancestral",
        "config": {
          "tree": zika.join("tree.nwk"),
          "alignment": [zika.join("aln.fasta.xz")],
        },
        "defer_start": true,
      })
      .to_string()
    }
  }
}
