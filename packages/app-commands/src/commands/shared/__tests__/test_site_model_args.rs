#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::argument_error;
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::partition::split::ModelSettings;

  const OPTIMIZE: [&str; 5] = ["optimize", "--tree", "t.nwk", "--alignment", "a.fasta"];
  const ANCESTRAL: [&str; 5] = ["ancestral", "--tree", "t.nwk", "--alignment", "a.fasta"];
  const TIMETREE: [&str; 7] = [
    "timetree",
    "--tree",
    "t.nwk",
    "--metadata",
    "m.tsv",
    "--alignment",
    "a.fasta",
  ];

  #[rustfmt::skip]
  #[rstest]
  #[case::no_rate_variation(&[],                          ModelSettings { gamma_categories: 1 })]
  #[case::four_categories(  &["--gamma-categories", "4"], ModelSettings { gamma_categories: 4 })]
  #[case::two_categories(   &["--gamma-categories", "2"], ModelSettings { gamma_categories: 2 })]
  #[trace]
  fn test_optimize_gamma_categories(#[case] args: &[&str], #[case] expected: ModelSettings) -> Result<(), Report> {
    let raw = TreetimeOptimizeArgsRaw::try_parse_from([&OPTIMIZE[..], args].concat())?;
    let actual = TreetimeOptimizeArgs::try_from(raw)?.sequences.models;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one(  "1")]
  #[case::zero( "0")]
  #[trace]
  fn test_too_few_gamma_categories_is_an_argument_error(#[case] categories: &str) {
    let raw = TreetimeOptimizeArgsRaw::try_parse_from([&OPTIMIZE[..], &["--gamma-categories", categories]].concat()).unwrap();
    let expected = (
      format!(
        "--gamma-categories needs 2 or more categories, but got {categories}: leave the setting out for a model \
         without rate variation"
      ),
      vec!["gamma_categories".to_owned()],
    );
    assert_eq!(Some(expected), argument_error(TreetimeOptimizeArgs::try_from(raw)));
  }

  #[test]
  fn test_ancestral_gamma_categories_reach_the_sequence_inputs() -> Result<(), Report> {
    let raw = TreetimeAncestralArgsRaw::try_parse_from([&ANCESTRAL[..], &["--gamma-categories", "4"]].concat())?;
    let actual = TreetimeAncestralArgs::try_from(raw)?.sequences.models;
    assert_eq!(ModelSettings { gamma_categories: 4 }, actual);
    Ok(())
  }

  #[test]
  fn test_ancestral_gamma_categories_with_parsimony_is_an_argument_error() {
    let raw = TreetimeAncestralArgsRaw::try_parse_from(
      [
        &ANCESTRAL[..],
        &["--method-anc", "parsimony", "--gamma-categories", "4"],
      ]
      .concat(),
    )
    .unwrap();
    let expected = (
      "--gamma-categories applies only with --method-anc marginal: parsimony reconstruction uses no substitution \
       model and no site rates"
        .to_owned(),
      vec!["gamma_categories".to_owned(), "method_anc".to_owned()],
    );
    assert_eq!(Some(expected), argument_error(TreetimeAncestralArgs::try_from(raw)));
  }

  #[test]
  fn test_timetree_gamma_categories_reach_the_sequence_inputs() -> Result<(), Report> {
    let raw = TreetimeTimetreeArgsRaw::try_parse_from([&TIMETREE[..], &["--gamma-categories", "4"]].concat())?;
    let actual = TreetimeTimetreeArgs::try_from(raw)?.sequences.models;
    assert_eq!(ModelSettings { gamma_categories: 4 }, actual);
    Ok(())
  }

  #[test]
  fn test_timetree_gamma_categories_with_input_branch_lengths_is_an_argument_error() {
    let raw = TreetimeTimetreeArgsRaw::try_parse_from(
      [
        &TIMETREE[..],
        &["--branch-length-mode", "input", "--gamma-categories", "4"],
      ]
      .concat(),
    )
    .unwrap();
    let expected = (
      "--gamma-categories requires re-estimated branch lengths (--branch-length-mode marginal), but \
       --branch-length-mode is input"
        .to_owned(),
      vec!["gamma_categories".to_owned(), "branch_length_mode".to_owned()],
    );
    assert_eq!(Some(expected), argument_error(TreetimeTimetreeArgs::try_from(raw)));
  }
}
