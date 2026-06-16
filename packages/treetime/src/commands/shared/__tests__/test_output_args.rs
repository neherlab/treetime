#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgs;
  use crate::commands::clock::args::TreetimeClockArgs;
  use crate::commands::mugration::args::TreetimeMugrationArgs;
  use crate::commands::optimize::args::TreetimeOptimizeArgs;
  use crate::commands::prune::args::TreetimePruneArgs;
  use crate::commands::shared::output::{
    AncestralOutputSelection, ClockOutputSelection, MugrationOutputSelection, NwkStyleArg, TimetreeOutputSelection,
  };
  use crate::commands::timetree::args::TreetimeTimetreeArgs;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use std::path::Path;

  #[test]
  fn test_output_args_output_all_parses() {
    let args =
      TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-all=/tmp/out"]).unwrap();
    assert_eq!(args.output.output_all.as_deref(), Some(Path::new("/tmp/out")));
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
      args.output_selection,
      vec![
        AncestralOutputSelection::Nwk,
        AncestralOutputSelection::Nexus,
        AncestralOutputSelection::GraphJson
      ]
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
    assert_eq!(args.output.output_tree_nwk.as_deref(), Some(Path::new("/tmp/my.nwk")));
  }

  #[test]
  fn test_output_args_per_file_tree_nexus_parses() {
    let args =
      TreetimeOptimizeArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-tree-nexus=/tmp/tree.nexus"])
        .unwrap();
    assert_eq!(
      args.output.output_tree_nexus.as_deref(),
      Some(Path::new("/tmp/tree.nexus"))
    );
  }

  #[test]
  fn test_output_args_nwk_style_parses_csv() {
    let args = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-tree-nwk=/tmp/my.nwk",
      "--output-nwk-style=plain,beast",
    ])
    .unwrap();
    assert_eq!(
      args.output.output_nwk_style,
      vec![NwkStyleArg::Plain, NwkStyleArg::Beast]
    );
  }

  #[test]
  fn test_output_args_short_flag_o_for_output_all() {
    let args = TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "-O", "/tmp/out"]).unwrap();
    assert_eq!(args.output.output_all.as_deref(), Some(Path::new("/tmp/out")));
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
    assert_eq!(args.output_selection, vec![AncestralOutputSelection::All]);
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
  fn test_output_args_rejects_out_of_command_selection_variant() {
    // `tracelog` is timetree-only; the ancestral selection enum has no such variant.
    let result = TreetimeAncestralArgs::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=tracelog",
    ]);
    assert!(
      result.is_err(),
      "timetree-only selection variant should be rejected on ancestral"
    );
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
      Some(Path::new("/tmp/node.json"))
    );
  }

  #[test]
  fn test_ancestral_output_gtr_parses() {
    let args =
      TreetimeAncestralArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-gtr=/tmp/gtr.json"]).unwrap();
    assert_eq!(args.output_gtr.as_deref(), Some(Path::new("/tmp/gtr.json")));
  }

  #[test]
  fn test_timetree_output_clock_model_parses() {
    let args = TreetimeTimetreeArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-clock-model=/tmp/clock.json",
    ])
    .unwrap();
    assert_eq!(args.output_clock_model.as_deref(), Some(Path::new("/tmp/clock.json")));
  }

  #[test]
  fn test_timetree_output_confidence_tsv_parses() {
    let args = TreetimeTimetreeArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-confidence-tsv=/tmp/conf.tsv",
    ])
    .unwrap();
    assert_eq!(args.output_confidence_tsv.as_deref(), Some(Path::new("/tmp/conf.tsv")));
  }

  #[test]
  fn test_timetree_output_tracelog_parses_with_alias() {
    let canonical =
      TreetimeTimetreeArgs::try_parse_from(["treetime", "--metadata=/dev/null", "--output-tracelog=/tmp/trace.csv"])
        .unwrap();
    assert_eq!(canonical.output_tracelog.as_deref(), Some(Path::new("/tmp/trace.csv")));

    let aliased =
      TreetimeTimetreeArgs::try_parse_from(["treetime", "--metadata=/dev/null", "--tracelog=/tmp/trace.csv"]).unwrap();
    assert_eq!(aliased.output_tracelog.as_deref(), Some(Path::new("/tmp/trace.csv")));
  }

  #[test]
  fn test_timetree_selection_confidence_tsv_parses() {
    let args = TreetimeTimetreeArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=confidence-tsv,tracelog",
    ])
    .unwrap();
    assert_eq!(
      args.output_selection,
      vec![
        TimetreeOutputSelection::ConfidenceTsv,
        TimetreeOutputSelection::Tracelog
      ]
    );
  }

  #[test]
  fn test_clock_output_clock_csv_parses() {
    let args =
      TreetimeClockArgs::try_parse_from(["treetime", "--metadata=/dev/null", "--output-clock-csv=/tmp/clock.csv"])
        .unwrap();
    assert_eq!(args.output_clock_csv.as_deref(), Some(Path::new("/tmp/clock.csv")));
  }

  #[test]
  fn test_clock_rejects_augur_node_data_selection() {
    // Clock does not produce augur node data; its selection enum has no such variant.
    let result = TreetimeClockArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=augur-node-data",
    ]);
    assert!(result.is_err(), "augur-node-data is not selectable on clock");
  }

  #[test]
  fn test_clock_selection_clock_model_parses() {
    let args = TreetimeClockArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=clock-model,clock-csv",
    ])
    .unwrap();
    assert_eq!(
      args.output_selection,
      vec![ClockOutputSelection::ClockModel, ClockOutputSelection::ClockCsv]
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
    assert_eq!(args.output_traits_csv.as_deref(), Some(Path::new("/tmp/traits.csv")));
  }

  #[test]
  fn test_mugration_output_confidence_csv_parses_with_alias() {
    let canonical = TreetimeMugrationArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--attribute=country",
      "--output-confidence-csv=/tmp/conf.csv",
    ])
    .unwrap();
    assert_eq!(
      canonical.output_confidence_csv.as_deref(),
      Some(Path::new("/tmp/conf.csv"))
    );

    let aliased = TreetimeMugrationArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--attribute=country",
      "--confidence=/tmp/conf.csv",
    ])
    .unwrap();
    assert_eq!(
      aliased.output_confidence_csv.as_deref(),
      Some(Path::new("/tmp/conf.csv"))
    );
  }

  #[test]
  fn test_mugration_selection_confidence_csv_parses() {
    let args = TreetimeMugrationArgs::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--attribute=country",
      "--output-all=/tmp/out",
      "--output-selection=confidence-csv",
    ])
    .unwrap();
    assert_eq!(args.output_selection, vec![MugrationOutputSelection::ConfidenceCsv]);
  }

  #[test]
  fn test_prune_output_gtr_parses() {
    let args =
      TreetimePruneArgs::try_parse_from(["treetime", "--tree=/dev/null", "--output-gtr=/tmp/gtr.json"]).unwrap();
    assert_eq!(args.output_gtr.as_deref(), Some(Path::new("/tmp/gtr.json")));
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
    assert_eq!(args.output.output_all.as_deref(), Some(Path::new("/tmp/out")));
    assert_eq!(args.output_selection, vec![AncestralOutputSelection::Nwk]);
    assert_eq!(
      args.output.output_tree_nexus.as_deref(),
      Some(Path::new("/tmp/custom.nexus"))
    );
    assert_eq!(
      args.output_augur_node_data.as_deref(),
      Some(Path::new("/tmp/custom-node.json"))
    );
  }
}
