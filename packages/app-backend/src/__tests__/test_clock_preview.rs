#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::sparse;
  use crate::clock_preview::{CLOCK_PREVIEW_TIP_LIMIT, ClockPreviewRequest, PreviewRoot, clock_preview};
  use crate::config::catalog::command_settings;
  use crate::results::auspice::run_auspice;
  use crate::results::run_results::{CommandResults, run_results};
  use crate::runs::manager::RunManager;
  use app_commands::command::AppCommand;
  use helpers::{clades, data, finished_run, star_dataset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  const NO_POSITIVE_RATE: &str = "The clock fit of the outlier filter requires a positive clock rate, also with \
    --allow-negative-rate. Run with --clock-filter=0 to skip the outlier filter: Clock rate is negative for all root \
    positions. The data may lack temporal signal. Please specify --clock-rate explicitly.";

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

  #[rustfmt::skip]
  #[rstest]
  #[case::identical_sampling_dates( "((a:0.1,b:0.2):0.05,(c:0.3,d:0.4):0.05);", Some("name\tdate\na\t2000\nb\t2000\nc\t2000\nd\t2000\n"), NO_POSITIVE_RATE)]
  #[case::unreadable_metadata(      "((a:0.1,b:0.2):0.05,(c:0.3,d:0.4):0.05);", None,                                                         "When reading dates: When opening file '<metadata>': No such file or directory (os error 2)")]
  #[case::no_positive_rate(         "(a:0.4,b:0.3,c:0.2,d:0.1);",               Some("name\tdate\na\t2000\nb\t2001\nc\t2002\nd\t2003\n"), NO_POSITIVE_RATE)]
  #[trace]
  fn test_clock_preview_draws_the_input_tree_when_the_regression_fails(
    #[case] newick: &str,
    #[case] metadata: Option<&str>,
    #[case] problem: &str,
  ) {
    let dir = tempdir().unwrap();
    let tree = dir.path().join("tree.nwk");
    let metadata_path = dir.path().join("metadata.tsv");
    fs::write(&tree, format!("{newick}\n")).unwrap();
    if let Some(metadata) = metadata {
      fs::write(&metadata_path, metadata).unwrap();
    }
    let problem = problem.replace("<metadata>", &metadata_path.to_string_lossy());
    let preview = clock_preview(&ClockPreviewRequest {
      command: AppCommand::Clock,
      config: sparse(json!({ "tree": tree, "metadata": metadata_path })),
      beyond_limit: false,
    })
    .unwrap();
    let tree = preview.tree.unwrap();

    assert_eq!(
      (4, true, PreviewRoot::Input, 4, Some(problem.as_str())),
      (
        preview.tips,
        preview.clock.is_none(),
        tree.root,
        clades(&tree.auspice.0.tree).iter().map(Vec::len).max().unwrap_or(0),
        preview.problem.as_deref()
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::clock_rate(   "clock_rate",   false)]
  #[case::clock_filter( "clock_filter", true)]
  #[case::reroot(       "reroot",       true)]
  #[case::tree(         "tree",         true)]
  #[case::coalescent(   "coalescent",   false)]
  #[case::confidence(   "confidence",   false)]
  #[case::output(       "output_all",   false)]
  #[trace]
  fn test_clock_preview_settings_of_the_time_tree(#[case] key: &str, #[case] read: bool) {
    let settings = command_settings(AppCommand::Timetree).unwrap().preview_settings;
    assert_eq!(read, settings.iter().any(|setting| setting == key));
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
