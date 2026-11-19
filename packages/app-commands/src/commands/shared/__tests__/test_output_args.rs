#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
  use crate::commands::clock::args::TreetimeClockArgsRaw;
  use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
  use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
  use crate::commands::prune::args::TreetimePruneArgsRaw;
  use crate::commands::shared::output_args::{
    AncestralOutputSelection, ClockOutputSelection, MugrationOutputSelection, NewickAnnotationsArg,
    TimetreeOutputSelection,
  };
  use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
  use clap::Parser;
  use clap::error::ErrorKind;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::Path;

  #[test]
  fn test_output_args_output_all_parses() {
    let args =
      TreetimeAncestralArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-all=/tmp/out"]).unwrap();
    assert_eq!(Some(Path::new("/tmp/out")), args.output.output_all.as_deref());
  }

  #[test]
  fn test_output_args_output_selection_parses_csv() {
    let args = TreetimeAncestralArgsRaw::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=nwk,nexus,graph-json",
    ])
    .unwrap();
    assert_eq!(
      args.output_selection.output_selection,
      vec![
        AncestralOutputSelection::Nwk,
        AncestralOutputSelection::Nexus,
        AncestralOutputSelection::GraphJson
      ]
    );
  }

  #[test]
  fn test_output_args_output_selection_requires_output_all() {
    let result = TreetimeAncestralArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-selection=nwk"]);
    assert_eq!(
      Err(ErrorKind::MissingRequiredArgument),
      result.map(|_| ()).map_err(|error| error.kind())
    );
  }

  #[test]
  fn test_output_args_per_file_tree_nwk_parses() {
    let args =
      TreetimeAncestralArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-tree-nwk=/tmp/my.nwk"])
        .unwrap();
    assert_eq!(Some(Path::new("/tmp/my.nwk")), args.output.output_tree_nwk.as_deref());
  }

  #[test]
  fn test_output_args_per_file_tree_nexus_parses() {
    let args =
      TreetimeOptimizeArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-tree-nexus=/tmp/tree.nexus"])
        .unwrap();
    assert_eq!(
      Some(Path::new("/tmp/tree.nexus")),
      args.output.output_tree_nexus.as_deref()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default(&[],                              NewickAnnotationsArg::Plain)]
  #[case::beast(  &["--newick-annotations=beast"],  NewickAnnotationsArg::Beast)]
  #[case::nhx(    &["--newick-annotations", "nhx"], NewickAnnotationsArg::Nhx)]
  #[trace]
  fn test_output_args_newick_annotations_parse(#[case] flags: &[&str], #[case] expected: NewickAnnotationsArg) {
    let args = TreetimeAncestralArgsRaw::try_parse_from(
      ["treetime", "--tree=/dev/null", "--output-tree-nwk=/tmp/my.nwk"].iter().chain(flags),
    )
    .unwrap();
    assert_eq!(expected, args.output.newick_annotations);
  }

  #[test]
  fn test_output_args_newick_annotations_take_one_value() {
    let result = TreetimeAncestralArgsRaw::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-tree-nwk=/tmp/my.nwk",
      "--newick-annotations=plain,beast",
    ]);
    assert_eq!(
      Err(ErrorKind::InvalidValue),
      result.map(|_| ()).map_err(|error| error.kind())
    );
  }

  #[test]
  fn test_output_args_short_flag_o_for_output_all() {
    let args = TreetimeAncestralArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "-O", "/tmp/out"]).unwrap();
    assert_eq!(Some(Path::new("/tmp/out")), args.output.output_all.as_deref());
  }

  #[test]
  fn test_output_args_output_selection_all_variant_parses() {
    let args = TreetimeAncestralArgsRaw::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=all",
    ])
    .unwrap();
    assert_eq!(
      vec![AncestralOutputSelection::All],
      args.output_selection.output_selection
    );
  }

  #[test]
  fn test_output_args_rejects_out_of_command_selection_variant() {
    let result = TreetimeAncestralArgsRaw::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=tracelog",
    ]);
    assert_eq!(
      Err(ErrorKind::InvalidValue),
      result.map(|_| ()).map_err(|error| error.kind())
    );
  }

  #[test]
  fn test_ancestral_output_node_data_parses() {
    let args =
      TreetimeAncestralArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-node-data=/tmp/node.json"])
        .unwrap();
    assert_eq!(Some(Path::new("/tmp/node.json")), args.output_node_data.as_deref());
  }

  #[test]
  fn test_ancestral_output_model_parses() {
    let args =
      TreetimeAncestralArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-model=/tmp/gtr.json"])
        .unwrap();
    assert_eq!(
      Some(Path::new("/tmp/gtr.json")),
      args.gtr_output.output_model.as_deref()
    );
  }

  #[test]
  fn test_timetree_output_clock_model_parses() {
    let args = TreetimeTimetreeArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-clock-model=/tmp/clock.json",
    ])
    .unwrap();
    assert_eq!(
      Some(Path::new("/tmp/clock.json")),
      args.clock_output.output_clock_model.as_deref()
    );
  }

  #[test]
  fn test_timetree_output_dates_parses() {
    let args =
      TreetimeTimetreeArgsRaw::try_parse_from(["treetime", "--metadata=/dev/null", "--output-dates=/tmp/dates.tsv"])
        .unwrap();
    assert_eq!(Some(Path::new("/tmp/dates.tsv")), args.output_dates.as_deref());
  }

  #[test]
  fn test_timetree_output_tracelog_parses() {
    let args =
      TreetimeTimetreeArgsRaw::try_parse_from(["treetime", "--metadata=/dev/null", "--output-tracelog=/tmp/trace.csv"])
        .unwrap();
    assert_eq!(Some(Path::new("/tmp/trace.csv")), args.output_tracelog.as_deref());
  }

  #[test]
  fn test_timetree_rejects_tracelog_without_output_prefix() {
    let result =
      TreetimeTimetreeArgsRaw::try_parse_from(["treetime", "--metadata=/dev/null", "--tracelog=/tmp/trace.csv"]);
    assert_eq!(
      Err(ErrorKind::UnknownArgument),
      result.map(|_| ()).map_err(|error| error.kind())
    );
  }

  #[test]
  fn test_timetree_selection_dates_parses() {
    let args = TreetimeTimetreeArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=dates,tracelog",
    ])
    .unwrap();
    assert_eq!(
      vec![TimetreeOutputSelection::Dates, TimetreeOutputSelection::Tracelog],
      args.output_selection.output_selection
    );
  }

  #[test]
  fn test_clock_output_clock_csv_parses() {
    let args =
      TreetimeClockArgsRaw::try_parse_from(["treetime", "--metadata=/dev/null", "--output-rtt=/tmp/rtt.tsv"]).unwrap();
    assert_eq!(Some(Path::new("/tmp/rtt.tsv")), args.clock_output.output_rtt.as_deref());
  }

  #[test]
  fn test_clock_rejects_node_data_selection() {
    let result = TreetimeClockArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=node-data",
    ]);
    assert_eq!(
      Err(ErrorKind::InvalidValue),
      result.map(|_| ()).map_err(|error| error.kind())
    );
  }

  #[test]
  fn test_clock_selection_clock_model_parses() {
    let args = TreetimeClockArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=clock-model,rtt",
    ])
    .unwrap();
    assert_eq!(
      vec![ClockOutputSelection::ClockModel, ClockOutputSelection::Rtt],
      args.output_selection.output_selection
    );
  }

  #[test]
  fn test_clock_output_chart_flags_parse() {
    let args = TreetimeClockArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--output-rtt-chart=/tmp/chart.png",
      "--output-all=/tmp/out",
      "--output-selection=rtt-chart",
    ])
    .unwrap();
    assert_eq!(
      (args.output_rtt_chart.as_deref(), args.output_selection.output_selection),
      (Some(Path::new("/tmp/chart.png")), vec![ClockOutputSelection::RttChart])
    );
  }

  #[test]
  fn test_mugration_output_traits_csv_parses() {
    let args = TreetimeMugrationArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--attribute=country",
      "--output-traits=/tmp/traits.csv",
    ])
    .unwrap();
    assert_eq!(Some(Path::new("/tmp/traits.csv")), args.output_traits.as_deref());
  }

  #[test]
  fn test_mugration_selection_confidence_csv_parses() {
    let args = TreetimeMugrationArgsRaw::try_parse_from([
      "treetime",
      "--metadata=/dev/null",
      "--attribute=country",
      "--output-all=/tmp/out",
      "--output-selection=trait-probabilities",
    ])
    .unwrap();
    assert_eq!(
      vec![MugrationOutputSelection::TraitProbabilities],
      args.output_selection.output_selection
    );
  }

  #[test]
  fn test_prune_rejects_the_model_output() {
    let result =
      TreetimePruneArgsRaw::try_parse_from(["treetime", "--tree=/dev/null", "--output-model=/tmp/model.json"]);
    assert_eq!(
      Err(ErrorKind::UnknownArgument),
      result.map(|_| ()).map_err(|error| error.kind())
    );
  }

  #[test]
  fn test_combined_output_all_and_per_file_parse() {
    let args = TreetimeAncestralArgsRaw::try_parse_from([
      "treetime",
      "--tree=/dev/null",
      "--output-all=/tmp/out",
      "--output-selection=nwk",
      "--output-tree-nexus=/tmp/custom.nexus",
      "--output-node-data=/tmp/custom-node.json",
    ])
    .unwrap();
    assert_eq!(Some(Path::new("/tmp/out")), args.output.output_all.as_deref());
    assert_eq!(
      vec![AncestralOutputSelection::Nwk],
      args.output_selection.output_selection
    );
    assert_eq!(
      Some(Path::new("/tmp/custom.nexus")),
      args.output.output_tree_nexus.as_deref()
    );
    assert_eq!(
      Some(Path::new("/tmp/custom-node.json")),
      args.output_node_data.as_deref()
    );
  }
}
