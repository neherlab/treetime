#[cfg(test)]
mod tests {
  use crate::runs::downloads::{download_name, output_label};
  use crate::runs::manager::RunManager;
  use app_commands::command::AppCommand;
  use app_output::output_plan::OutputSelection;
  use helpers::titled_run;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;

  #[test]
  fn test_downloads_of_two_runs_with_the_same_title_have_different_names() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let first = titled_run(&runs, "Ebola/100, time tree: skyline");
    let second = titled_run(&runs, "Ebola/100, time tree: skyline");

    let names = [&first, &second].map(|record| download_name(record, "timetree.nwk"));

    let short_id = |id: &str| id.chars().take(8).collect::<String>();
    assert_eq!(
      [
        format!(
          "ebola-100-time-tree-skyline-{}.timetree.nwk",
          short_id(first.id.as_str())
        ),
        format!(
          "ebola-100-time-tree-skyline-{}.timetree.nwk",
          short_id(second.id.as_str())
        ),
      ],
      names
    );
    assert_ne!(names[0], names[1]);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::symbols_only(   "*** / ***",             "treetime-<id>.partition-tree.nwk")]
  #[case::accented(       "Zürich Ñandú",          "zürich-ñandú-<id>.partition-tree.nwk")]
  #[case::long_title(     &"word ".repeat(30),     &format!("{}-<id>.partition-tree.nwk", "word-".repeat(12).trim_end_matches('-')))]
  #[trace]
  fn test_downloads_name_keeps_letters_and_digits_of_the_title(#[case] title: &str, #[case] expected: &str) {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = titled_run(&runs, title);

    let id = record.id.as_str().chars().take(8).collect::<String>();
    assert_eq!(expected.replace("<id>", &id), download_name(&record, "partition/tree.nwk"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::titled(      "Ebola/100, time tree: skyline", "ebola-100-time-tree-skyline-<id>.zip")]
  #[case::symbols_only("*** / ***",                     "treetime-<id>.zip")]
  #[trace]
  fn test_downloads_package_name_has_the_stem_of_the_output_names(#[case] title: &str, #[case] expected: &str) {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = titled_run(&runs, title);

    let id = record.id.as_str().chars().take(8).collect::<String>();
    assert_eq!(expected.replace("<id>", &id), record.summary().package_name);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dated_tree(     AppCommand::Timetree, Some(OutputSelection::Nwk),      "timetree.nwk",          "Dated tree (Newick, branch lengths in years)")]
  #[case::rerooted_tree(  AppCommand::Clock,    Some(OutputSelection::Nwk),      "clock.nwk",             "Rerooted tree (Newick)")]
  #[case::ancestral_fasta(AppCommand::Ancestral,Some(OutputSelection::Sequences),"ancestral.1.sequences.fasta", "Ancestral sequences (FASTA)")]
  #[case::png_chart(      AppCommand::Clock,    Some(OutputSelection::RttChart), "clock.rtt.png",         "Root-to-tip chart (PNG)")]
  #[case::unplanned_file( AppCommand::Timetree, None,                            "logs/extra.txt",        "extra.txt")]
  #[trace]
  fn test_downloads_label_names_the_content_and_format(
    #[case] command: AppCommand,
    #[case] kind: Option<OutputSelection>,
    #[case] path: &str,
    #[case] expected: &str,
  ) {
    assert_eq!(expected, output_label(command, kind, path));
  }

  mod helpers {
    use crate::__tests__::test_support::tests::sparse;
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunRecord, UpdateRunRequest};
    use app_commands::command::AppCommand;
    use serde_json::json;
    use std::path::Path;

    pub(super) fn titled_run(runs: &RunManager, title: &str) -> RunRecord {
      let created = runs
        .create(&CreateRunRequest {
          command: AppCommand::Clock,
          config: sparse(json!({ "tree": Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20/tree.nwk") })),
          defer_start: true,
          source_run_id: None,
          example: None,
        })
        .unwrap();
      runs
        .update(
          &created.id,
          UpdateRunRequest {
            title: Some(title.to_owned()),
            pinned: None,
          },
        )
        .unwrap();
      runs.get(&created.id).unwrap()
    }
  }
}
