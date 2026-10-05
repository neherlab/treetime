#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::path::Path;

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree_defaults( AppCommand::Timetree,  json!(null),             json!(["nwk", "nexus", "auspice", "augur-node-data", "gtr", "reconstructed-nuc-fasta", "clock-model", "coalescent-tsv", "tracelog", "clock-csv"]))]
  #[case::clock_defaults(    AppCommand::Clock,     json!(null),             json!(["nwk", "nexus", "clock-model", "clock-csv", "clock-chart-svg", "clock-chart-png", "auspice"]))]
  #[case::ancestral_defaults(AppCommand::Ancestral, json!(null),             json!(["nwk", "nexus", "augur-node-data", "gtr", "reconstructed-nuc-fasta", "auspice"]))]
  #[case::mugration_defaults(AppCommand::Mugration, json!(null),             json!(["nwk", "nexus", "augur-node-data", "gtr", "traits-csv", "auspice"]))]
  #[case::optimize_defaults( AppCommand::Optimize,  json!(null),             json!(["nwk", "nexus", "augur-node-data", "gtr", "auspice"]))]
  #[case::prune_defaults(    AppCommand::Prune,     json!(null),             json!(["nwk", "nexus", "gtr", "auspice"]))]
  #[case::prune_chosen(      AppCommand::Prune,     json!(["mat-pb"]),        json!(["mat-pb", "auspice"]))]
  #[case::timetree_chosen(   AppCommand::Timetree,  json!(["nwk"]),          json!(["nwk", "auspice", "tracelog", "clock-csv"]))]
  #[case::timetree_has_one(  AppCommand::Timetree,  json!(["tracelog"]),     json!(["tracelog", "auspice", "clock-csv"]))]
  #[case::clock_all(         AppCommand::Clock,     json!(["all"]),          json!(["all"]))]
  #[trace]
  fn test_prepare_run_adds_the_outputs_the_views_need(
    #[case] command: AppCommand,
    #[case] chosen: Value,
    #[case] expected: Value,
  ) {
    let mut config = match command {
      AppCommand::Timetree | AppCommand::Clock => json!({ "tree": "t.nwk", "metadata": "m.tsv" }),
      AppCommand::Mugration => json!({ "tree": "t.nwk", "metadata": "m.tsv", "attribute": "country" }),
      AppCommand::Ancestral | AppCommand::Optimize | AppCommand::Prune => json!({ "tree": "t.nwk" }),
    };
    if !chosen.is_null() {
      config["output_selection"] = chosen;
    }
    let prepared = command.prepare_run(&config, Path::new("/runs/r/out")).unwrap();
    assert_eq!(
      (expected, json!("/runs/r/out")),
      (prepared.config["output_selection"].clone(), prepared.config["output_all"].clone())
    );
  }

  #[test]
  fn test_prepare_run_drops_client_output_paths() {
    let config = json!({
      "tree": "t.nwk",
      "metadata": "m.tsv",
      "output_all": "/etc",
      "output_tree_nwk": "/etc/passwd",
      "output_clock_model": "../x.json",
      "output_nwk_style": ["beast"],
    });
    let prepared = AppCommand::Clock
      .prepare_run(&config, Path::new("/runs/r/out"))
      .unwrap();
    assert_eq!(
      (Some(json!("/runs/r/out")), None, None, Some(json!(["beast"]))),
      (
        prepared.config.get("output_all").cloned(),
        prepared.config.get("output_tree_nwk").cloned(),
        prepared.config.get("output_clock_model").cloned(),
        prepared.config.get("output_nwk_style").cloned()
      )
    );
  }

  #[test]
  fn test_prepare_run_lists_settings_changed_from_the_defaults() {
    let config = json!({
      "tree": "t.nwk",
      "metadata": "m.tsv",
      "clock_filter": 2.0,
      "keep_root": true,
      "branch_split": { "n_points": 7 },
    });
    let prepared = AppCommand::Clock
      .prepare_run(&config, Path::new("/runs/r/out"))
      .unwrap();
    let mut changed = prepared.changed_settings;
    changed.sort();
    assert_eq!(vec!["branch_split.n_points", "clock_filter", "keep_root"], changed);
  }
}
