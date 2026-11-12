#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::sparse;
  use crate::clock_preview::{CLOCK_PREVIEW_TIP_LIMIT, ClockPreviewRequest, PreviewRoot, clock_preview};
  use crate::results::auspice::run_auspice;
  use crate::results::run_results::{CommandResults, run_results};
  use crate::runs::manager::RunManager;
  use app_commands::command::AppCommand;
  use helpers::{clades, data, finished_run, star_dataset};
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

    let run_tree = run_auspice(&runs, &id).unwrap();
    let tree = preview.tree.unwrap();

    assert_eq!(
      (
        20,
        CLOCK_PREVIEW_TIP_LIMIT,
        Some(expected),
        PreviewRoot::Regression,
        clades(&run_tree.0.tree)
      ),
      (
        preview.tips,
        preview.tip_limit,
        preview.clock,
        tree.root,
        clades(&tree.auspice.0.tree)
      )
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
      (tips, CLOCK_PREVIEW_TIP_LIMIT, computed, computed),
      (
        preview.tips,
        preview.tip_limit,
        preview.clock.is_some(),
        preview.tree.is_some()
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_metadata( AppCommand::Timetree,  json!({ "tree": data("zika/20/tree.nwk") }))]
  #[case::without_dates(    AppCommand::Ancestral, json!({ "tree": data("zika/20/tree.nwk"), "alignment": [data("zika/20/aln.fasta.xz")] }))]
  #[trace]
  fn test_clock_preview_without_a_regression_draws_the_input_tree(#[case] command: AppCommand, #[case] draft: Value) {
    let preview = clock_preview(&ClockPreviewRequest {
      command,
      config: sparse(draft),
      beyond_limit: false,
    })
    .unwrap();
    let tree = preview.tree.unwrap();
    let root = &tree.auspice.0.tree;

    assert_eq!(
      (20, true, PreviewRoot::Input, 20, Some(0.0)),
      (
        preview.tips,
        preview.clock.is_none(),
        tree.root,
        clades(root).iter().map(Vec::len).max().unwrap_or(0),
        root.node_attrs.div
      )
    );
  }

  #[test]
  fn test_clock_preview_needs_a_tree() {
    assert_error!(
      clock_preview(&ClockPreviewRequest {
        command: AppCommand::Timetree,
        config: sparse(json!({ "metadata": data("zika/20/metadata.tsv") })),
        beyond_limit: false,
      }),
      "a clock preview needs a tree"
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
    use std::collections::BTreeSet;
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use treetime_io::auspice_types::AuspiceTreeNode;
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

    pub(super) fn clades(node: &AuspiceTreeNode) -> BTreeSet<Vec<String>> {
      let mut clades = BTreeSet::new();
      collect_clades(node, &mut clades);
      clades
    }

    fn collect_clades(node: &AuspiceTreeNode, clades: &mut BTreeSet<Vec<String>>) -> Vec<String> {
      if node.children.is_empty() {
        return vec![node.name.clone()];
      }
      let tips = node
        .children
        .iter()
        .flat_map(|child| collect_clades(child, clades))
        .sorted()
        .collect_vec();
      clades.insert(tips.clone());
      tips
    }

    pub(super) fn data(path: &str) -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(path)
    }
  }
}
