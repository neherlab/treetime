#[cfg(test)]
mod tests {
  use crate::runs::manager::RunManager;
  use crate::runs::record::UpdateRunRequest;
  use helpers::{prepared, zika};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use tempfile::tempdir;

  #[test]
  fn test_run_name_of_a_default_run_names_the_example_and_the_analysis() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let run = prepared(&runs, zika(&json!({})), Some("zika/20"), None);
    assert_eq!("zika/20, time tree", run.title);
  }

  #[test]
  fn test_run_name_of_own_data_names_the_tree_file() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let run = prepared(&runs, zika(&json!({})), None, None);
    assert_eq!("tree.nwk, time tree", run.title);
  }

  #[test]
  fn test_run_name_of_a_variant_names_its_change_against_the_source() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let source = prepared(&runs, zika(&json!({ "confidence": true })), Some("zika/20"), None);
    let variant = prepared(
      &runs,
      zika(&json!({ "confidence": true, "coalescent": "skyline" })),
      Some("zika/20"),
      Some(&source.id),
    );
    assert_eq!(
      ("zika/20, time tree: confidence", "zika/20, time tree: skyline"),
      (source.title.as_str(), variant.title.as_str())
    );
  }

  #[test]
  fn test_run_name_a_variant_with_the_inputs_of_its_source_keeps_the_example() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let source = prepared(&runs, zika(&json!({})), Some("zika/20"), None);
    let variant = prepared(&runs, zika(&json!({ "confidence": true })), None, Some(&source.id));
    assert_eq!(
      (Some("zika/20"), "zika/20, time tree: confidence", Some("zika/20")),
      (
        variant.example.as_deref(),
        variant.title.as_str(),
        variant.dataset.as_deref()
      )
    );
  }

  #[test]
  fn test_run_name_a_variant_with_other_inputs_names_its_tree_file() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let source = prepared(&runs, zika(&json!({})), Some("zika/20"), None);
    let mut config = zika(&json!({}));
    config["alignment"] = json!([]);
    config["sequence_length"] = json!(10_000);
    let variant = prepared(&runs, config, None, Some(&source.id));
    assert_eq!(
      (None, "tree.nwk"),
      (variant.example.as_deref(), variant.summary().dataset.unwrap().as_str())
    );
  }

  #[test]
  fn test_run_name_a_run_is_prepared_after_its_start_only() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = helpers::create(&runs, zika(&json!({})), None, None);
    let started = helpers::start(&runs, &created.id);
    assert_eq!((false, true), (created.prepared, started.prepared));
  }

  #[test]
  fn test_run_name_names_three_changes_and_counts_the_others() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let run = prepared(
      &runs,
      zika(&json!({
        "clock_filter": 3.0,
        "max_iter": 4,
        "confidence": true,
        "keep_polytomies": true,
        "coalescent": "opt",
      })),
      Some("zika/20"),
      None,
    );
    assert_eq!(
      "zika/20, time tree: clock filter 3, iterations 4, confidence +2",
      run.title
    );
  }

  #[test]
  fn test_run_name_leaves_out_settings_that_are_no_assumptions() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let run = prepared(
      &runs,
      zika(&json!({ "seed": 3, "date_format": "%Y", "output_selection": ["nwk"], "newick_annotations": "nhx" })),
      Some("zika/20"),
      None,
    );
    assert_eq!("zika/20, time tree", run.title);
  }

  #[test]
  fn test_run_name_of_a_switched_off_setting_says_off() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let source = prepared(&runs, zika(&json!({ "confidence": true })), None, None);
    let variant = prepared(&runs, zika(&json!({})), None, Some(&source.id));
    assert_eq!("tree.nwk, time tree: confidence off", variant.title);
  }

  #[test]
  fn test_run_name_of_a_renamed_run_stays_when_the_run_starts() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = helpers::create(&runs, zika(&json!({ "confidence": true })), Some("zika/20"), None);
    runs
      .update(
        &created.id,
        UpdateRunRequest {
          title: Some("my run".to_owned()),
          pinned: None,
        },
      )
      .unwrap();
    let run = helpers::start(&runs, &created.id);
    assert_eq!(("my run", true), (run.title.as_str(), run.renamed));
  }

  mod helpers {
    use crate::__tests__::test_support::tests::{sparse, unconfined};
    use crate::job::JobId;
    use crate::runs::lifecycle::Prepared;
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunRecord, StartRunRequest};
    use app_commands::command::AppCommand;
    use serde_json::{Value, json};
    use std::path::Path;
    use std::sync::Arc;

    pub(super) fn zika(settings: &Value) -> Value {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      let mut config = json!({
        "tree": zika.join("tree.nwk"),
        "metadata": zika.join("metadata.tsv"),
        "alignment": [zika.join("aln.fasta.xz")],
      });
      if let (Value::Object(config), Value::Object(settings)) = (&mut config, settings) {
        config.extend(settings.clone());
      }
      config
    }

    pub(super) fn prepared(
      runs: &Arc<RunManager>,
      config: Value,
      example: Option<&str>,
      source: Option<&JobId>,
    ) -> RunRecord {
      let created = create(runs, config, example, source);
      start(runs, &created.id)
    }

    pub(super) fn create(
      runs: &Arc<RunManager>,
      config: Value,
      example: Option<&str>,
      source: Option<&JobId>,
    ) -> RunRecord {
      runs
        .create(&CreateRunRequest {
          command: AppCommand::Timetree,
          config: sparse(config),
          defer_start: true,
          source_run_id: source.cloned(),
          example: example.map(ToOwned::to_owned),
        })
        .unwrap()
    }

    pub(super) fn start(runs: &Arc<RunManager>, id: &JobId) -> RunRecord {
      let started = runs.start(id, &StartRunRequest::default(), unconfined()).unwrap();
      assert!(matches!(started.prepare(), Prepared::Ready(_)), "the run prepares");
      runs.get(id).unwrap()
    }
  }
}
