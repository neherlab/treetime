#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgs;
  use crate::commands::clock::args::TreetimeClockArgs;
  use crate::commands::mugration::args::TreetimeMugrationArgs;
  use crate::commands::optimize::args::TreetimeOptimizeArgs;
  use crate::commands::prune::args::TreetimePruneArgs;
  use crate::commands::shared::output::OutputSelection;
  use crate::commands::timetree::args::TreetimeTimetreeArgs;
  use clap::Parser;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_output_args_output_all_parses() {
    let args =
      TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-all=/tmp/out"]).unwrap();
    assert_eq!(
      args.output.output_all.as_deref(),
      Some(std::path::Path::new("/tmp/out"))
    );
  }

  #[test]
  fn test_output_args_output_selection_parses_csv() {
    let args = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=nwk,nexus,graph-json",
    ])
    .unwrap();
    assert_eq!(
      args.output.output_selection,
      vec![OutputSelection::Nwk, OutputSelection::Nexus, OutputSelection::GraphJson]
    );
  }

  #[test]
  fn test_output_args_output_selection_requires_output_all() {
    let result = TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-selection=nwk"]);
    assert!(result.is_err(), "--output-selection without --output-all should fail");
  }

  #[test]
  fn test_output_args_per_file_tree_nwk_parses() {
    let args =
      TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-tree-nwk=/tmp/my.nwk"]).unwrap();
    assert_eq!(
      args.output.output_tree_nwk.as_deref(),
      Some(std::path::Path::new("/tmp/my.nwk"))
    );
  }

  #[test]
  fn test_output_args_per_file_tree_nexus_annotated_parses() {
    let args = TreetimeOptimizeArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-tree-nexus-annotated=/tmp/tree.annotated.nexus",
    ])
    .unwrap();
    assert_eq!(
      args.output.output_tree_nexus_annotated.as_deref(),
      Some(std::path::Path::new("/tmp/tree.annotated.nexus"))
    );
  }

  #[test]
  fn test_output_args_short_flag_o_for_output_all() {
    let args = TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "-O", "/tmp/out"]).unwrap();
    assert_eq!(
      args.output.output_all.as_deref(),
      Some(std::path::Path::new("/tmp/out"))
    );
  }

  #[test]
  fn test_output_args_output_selection_all_variant_parses() {
    let args = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=all",
    ])
    .unwrap();
    assert_eq!(args.output.output_selection, vec![OutputSelection::All]);
  }

  #[test]
  fn test_output_args_rejects_unknown_selection_value() {
    let result = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=unknown-format",
    ]);
    assert!(result.is_err(), "unknown --output-selection value should be rejected");
  }

  #[test]
  fn test_ancestral_output_augur_node_data_parses() {
    let args = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-augur-node-data=/tmp/node.json",
    ])
    .unwrap();
    assert_eq!(
      args.output_augur_node_data.as_deref(),
      Some(std::path::Path::new("/tmp/node.json"))
    );
  }

  #[test]
  fn test_ancestral_output_gtr_parses() {
    let args =
      TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-gtr=/tmp/gtr.json"]).unwrap();
    assert_eq!(args.output_gtr.as_deref(), Some(std::path::Path::new("/tmp/gtr.json")));
  }

  #[test]
  fn test_timetree_output_clock_model_parses() {
    let args = TreetimeTimetreeArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-clock-model=/tmp/clock.json",
    ])
    .unwrap();
    assert_eq!(
      args.output_clock_model.as_deref(),
      Some(std::path::Path::new("/tmp/clock.json"))
    );
  }

  #[test]
  fn test_clock_output_clock_csv_parses() {
    let args =
      TreetimeClockArgs::try_parse_from(["treetime", "--metadata=/dev/null", "--output-clock-csv=/tmp/clock.csv"])
        .unwrap();
    assert_eq!(
      args.output_clock_csv.as_deref(),
      Some(std::path::Path::new("/tmp/clock.csv"))
    );
  }

  #[test]
  fn test_mugration_output_traits_csv_parses() {
    let args = TreetimeMugrationArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--attribute=country",
      "--output-traits-csv=/tmp/traits.csv",
    ])
    .unwrap();
    assert_eq!(
      args.output_traits_csv.as_deref(),
      Some(std::path::Path::new("/tmp/traits.csv"))
    );
  }

  #[test]
  fn test_prune_output_gtr_parses() {
    let args =
      TreetimePruneArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-gtr=/tmp/gtr.json"]).unwrap();
    assert_eq!(args.output_gtr.as_deref(), Some(std::path::Path::new("/tmp/gtr.json")));
  }

  #[test]
  fn test_combined_output_all_and_per_file_parse() {
    let args = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=nwk",
      "--output-tree-nexus=/tmp/custom.nexus",
      "--output-augur-node-data=/tmp/custom-node.json",
    ])
    .unwrap();
    assert_eq!(
      args.output.output_all.as_deref(),
      Some(std::path::Path::new("/tmp/out"))
    );
    assert_eq!(args.output.output_selection, vec![OutputSelection::Nwk]);
    assert_eq!(
      args.output.output_tree_nexus.as_deref(),
      Some(std::path::Path::new("/tmp/custom.nexus"))
    );
    assert_eq!(
      args.output_augur_node_data.as_deref(),
      Some(std::path::Path::new("/tmp/custom-node.json"))
    );
  }
}
