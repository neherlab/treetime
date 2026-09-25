#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::path::Path;

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree_defaults( AppCommand::Timetree,  json!(null),             json!(["Nwk", "Nexus", "Auspice", "AugurNodeData", "Gtr", "ReconstructedNucFasta", "ClockModel", "CoalescentTsv", "Tracelog"]))]
  #[case::clock_defaults(    AppCommand::Clock,     json!(null),             json!(["Nwk", "Nexus", "ClockModel", "ClockCsv", "Auspice"]))]
  #[case::ancestral_defaults(AppCommand::Ancestral, json!(null),             json!(["Nwk", "Nexus", "AugurNodeData", "Gtr", "ReconstructedNucFasta", "Auspice"]))]
  #[case::mugration_defaults(AppCommand::Mugration, json!(null),             json!(["Nwk", "Nexus", "AugurNodeData", "Gtr", "TraitsCsv", "Auspice"]))]
  #[case::optimize_defaults( AppCommand::Optimize,  json!(null),             json!(["Nwk", "Nexus", "AugurNodeData", "Gtr", "Auspice"]))]
  #[case::prune_defaults(    AppCommand::Prune,     json!(null),             json!(["Nwk", "Nexus", "Gtr"]))]
  #[case::timetree_chosen(   AppCommand::Timetree,  json!(["Nwk"]),          json!(["Nwk", "Auspice", "Tracelog"]))]
  #[case::timetree_has_one(  AppCommand::Timetree,  json!(["Tracelog"]),     json!(["Tracelog", "Auspice"]))]
  #[case::clock_all(         AppCommand::Clock,     json!(["All"]),          json!(["All"]))]
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
      (json!("/runs/r/out"), json!(null), json!(null), json!(["beast"])),
      (
        prepared.config["output_all"].clone(),
        prepared.config["output_tree_nwk"].clone(),
        prepared.config["output_clock_model"].clone(),
        prepared.config["output_nwk_style"].clone()
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
