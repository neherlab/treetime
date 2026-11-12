#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::sparse;
  use crate::clock_preview::{CLOCK_PREVIEW_TIP_LIMIT, ClockPreview, ClockPreviewRequest, clock_preview};
  use crate::results::run_results::{CommandResults, run_results};
  use crate::runs::manager::RunManager;
  use app_commands::command::AppCommand;
  use helpers::{data, finished_run, star_dataset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::clock_draft(    AppCommand::Clock,    json!({ "tree": data("zika/20/tree.nwk"), "metadata": data("zika/20/metadata.tsv") }),
                                                json!({ "tree": data("zika/20/tree.nwk"), "metadata": data("zika/20/metadata.tsv") }))]
  #[case::timetree_draft( AppCommand::Timetree, json!({ "tree": data("ebola/20/tree.nwk"), "metadata": data("ebola/20/metadata.tsv"), "alignment": [data("ebola/20/aln.fasta.xz")], "clock_rate": 0.001, "clock_filter": 2.5, "coalescent": "skyline" }),
                                                json!({ "tree": data("ebola/20/tree.nwk"), "metadata": data("ebola/20/metadata.tsv"), "clock_filter": 2.5 }))]
  #[case::min_dev_root(   AppCommand::Timetree, json!({ "tree": data("ebola/20/tree.nwk"), "metadata": data("ebola/20/metadata.tsv"), "reroot": "min-dev" }),
                                                json!({ "tree": data("ebola/20/tree.nwk"), "metadata": data("ebola/20/metadata.tsv"), "reroot": "min-dev" }))]
  #[trace]
  fn test_clock_preview_equals_the_root_to_tip_data_of_a_clock_run(
    #[case] command: AppCommand,
    #[case] draft: Value,
    #[case] clock: Value,
  ) {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let id = finished_run(&runs, clock);
    let CommandResults::Clock(expected) = run_results(&runs, &id).unwrap().results else {
      panic!("a clock run has clock results");
    };

    let preview = clock_preview(&ClockPreviewRequest {
      command,
      config: sparse(draft),
      beyond_limit: false,
    })
    .unwrap();

    assert_eq!(
      ClockPreview {
        tips: 20,
        tip_limit: CLOCK_PREVIEW_TIP_LIMIT,
        clock: Some(expected),
      },
      preview
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::above_limit_on_request(true,  true)]
  #[case::above_limit(           false, false)]
  #[trace]
  fn test_clock_preview_of_a_tree_above_the_tip_limit_computes_only_on_request(
    #[case] beyond_limit: bool,
    #[case] computed: bool,
  ) {
    let dir = tempdir().unwrap();
    let tips = CLOCK_PREVIEW_TIP_LIMIT + 1;
    let (tree, metadata) = star_dataset(dir.path(), tips);
    let preview = clock_preview(&ClockPreviewRequest {
      command: AppCommand::Timetree,
      config: sparse(json!({ "tree": tree, "metadata": metadata })),
      beyond_limit,
    })
    .unwrap();
    assert_eq!(
      (tips, CLOCK_PREVIEW_TIP_LIMIT, computed),
      (preview.tips, preview.tip_limit, preview.clock.is_some())
    );
  }

  #[test]
  fn test_clock_preview_needs_the_metadata() {
    assert_error!(
      clock_preview(&ClockPreviewRequest {
        command: AppCommand::Timetree,
        config: sparse(json!({ "tree": data("zika/20/tree.nwk") })),
        beyond_limit: false,
      }),
      "a clock preview needs the metadata with the sampling dates"
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::sparse;
    use crate::job::JobId;
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunStatus, StartRunRequest};
    use app_commands::command::AppCommand;
    use itertools::Itertools;
    use pretty_assertions::assert_eq;
    use serde_json::Value;
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn finished_run(runs: &Arc<RunManager>, config: Value) -> JobId {
      let created = runs
        .create(&CreateRunRequest {
          command: AppCommand::Clock,
          config: sparse(config),
          defer_start: true,
          source_run_id: None,
          example: None,
        })
        .unwrap();
      runs
        .start(
          &created.id,
          &StartRunRequest::default(),
          Box::new(|_config: &mut Value| Ok(())),
        )
        .unwrap()
        .run();
      assert_eq!(RunStatus::Ok, runs.get(&created.id).unwrap().status);
      created.id
    }

    pub(super) fn star_dataset(dir: &Path, tips: usize) -> (PathBuf, PathBuf) {
      let tree = dir.join("tree.nwk");
      let metadata = dir.join("metadata.tsv");
      let years = (0..tips).map(|tip| tip % 50).collect_vec();
      let branches = years
        .iter()
        .enumerate()
        .map(|(tip, year)| format!("t{tip}:{}", 0.001 * (year + 1).to_f64_exact()))
        .join(",");
      let rows = years
        .iter()
        .enumerate()
        .map(|(tip, year)| format!("t{tip}\t{}\n", 2000 + year))
        .join("");
      fs::write(&tree, format!("({branches});\n")).unwrap();
      fs::write(&metadata, format!("name\tdate\n{rows}")).unwrap();
      (tree, metadata)
    }

    pub(super) fn data(path: &str) -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(path)
    }
  }
}
