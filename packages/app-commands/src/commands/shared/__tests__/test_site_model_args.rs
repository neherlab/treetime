#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::argument_error;
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use crate::commands::shared::partition_by::PartitionScheme;
  use crate::commands::shared::site_model::default_shared_params;
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::gtr::shared_params::SharedParams;
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
  #[case::no_rate_variation(&[],                          ModelSettings { gamma_categories: 1, shared: SharedParams::default() })]
  #[case::four_categories(  &["--gamma-categories", "4"], ModelSettings { gamma_categories: 4, shared: SharedParams::default() })]
  #[case::two_categories(   &["--gamma-categories", "2"], ModelSettings { gamma_categories: 2, shared: SharedParams::default() })]
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
    assert_eq!(
      ModelSettings {
        gamma_categories: 4,
        shared: SharedParams::default()
      },
      actual
    );
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
       model"
        .to_owned(),
      vec!["gamma_categories".to_owned(), "method_anc".to_owned()],
    );
    assert_eq!(Some(expected), argument_error(TreetimeAncestralArgs::try_from(raw)));
  }

  #[test]
  fn test_timetree_gamma_categories_reach_the_sequence_inputs() -> Result<(), Report> {
    let raw = TreetimeTimetreeArgsRaw::try_parse_from([&TIMETREE[..], &["--gamma-categories", "4"]].concat())?;
    let actual = TreetimeTimetreeArgs::try_from(raw)?.sequences.models;
    assert_eq!(
      ModelSettings {
        gamma_categories: 4,
        shared: SharedParams::default()
      },
      actual
    );
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

  #[rustfmt::skip]
  #[rstest]
  #[case::files(              None,                                 false, (false, false, false))]
  #[case::files_gamma(        None,                                 true,  (false, false, false))]
  #[case::codon_positions(    Some(PartitionScheme::CodonPositions), true,  (false, false, false))]
  #[case::genes(              Some(PartitionScheme::Genes),          false, (true,  true,  false))]
  #[case::genes_gamma(        Some(PartitionScheme::Genes),          true,  (true,  true,  true))]
  #[case::rate_categories(    Some(PartitionScheme::RateCategories), false, (true,  true,  false))]
  #[case::rate_categories_gamma(Some(PartitionScheme::RateCategories), true, (true,  true,  true))]
  #[trace]
  fn test_default_shared_params_follow_the_scheme(
    #[case] scheme: Option<PartitionScheme>,
    #[case] gamma: bool,
    #[case] (exchange_rates, frequencies, gamma_shape): (bool, bool, bool),
  ) {
    let expected = SharedParams {
      exchange_rates,
      frequencies,
      gamma_shape,
    };
    assert_eq!(expected, default_shared_params(scheme, gamma));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::genes_default(    &["--partition-by", "genes", "--annotation", "a.gff"],                                  (true,  true,  false))]
  #[case::genes_named_model(&["--partition-by", "genes", "--annotation", "a.gff", "--model", "jc69"],               (false, false, false))]
  #[case::genes_none(       &["--partition-by", "genes", "--annotation", "a.gff", "--shared-params", "none"],       (false, false, false))]
  #[case::codon_frequencies(&["--partition-by", "codon-positions", "--annotation", "a.gff", "--shared-params", "frequencies"], (false, true, false))]
  #[case::entries_explicit( &["--partition", "{sites: ['1-10']}", "--partition", "{}", "--shared-params", "exchange-rates,frequencies"], (true, true, false))]
  #[case::one_file_default( &[],                                                                                    (false, false, false))]
  #[trace]
  fn test_optimize_shared_params(
    #[case] args: &[&str],
    #[case] (exchange_rates, frequencies, gamma_shape): (bool, bool, bool),
  ) -> Result<(), Report> {
    let raw = TreetimeOptimizeArgsRaw::try_parse_from([&OPTIMIZE[..], args].concat())?;
    let expected = SharedParams {
      exchange_rates,
      frequencies,
      gamma_shape,
    };
    assert_eq!(expected, TreetimeOptimizeArgs::try_from(raw)?.sequences.models.shared);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::shape_without_gamma(&["--partition-by", "genes", "--annotation", "a.gff", "--shared-params", "gamma-shape"],            "--shared-params gamma-shape needs --gamma-categories: without rate categories, the partitions have no gamma shape", &["shared_params", "gamma_categories"])]
  #[case::none_with_others(   &["--partition-by", "genes", "--annotation", "a.gff", "--shared-params", "none,frequencies"],       "--shared-params takes 'none' alone, but got 'none', 'frequencies'", &["shared_params"])]
  #[case::none_one_file(      &["--shared-params", "none"],                                                                     "--shared-params applies only with two or more partitions whose model is infer, but the run has 1 partition: without --partition or --partition-by, each --alignment file is one partition", &["shared_params", "alignment"])]
  #[case::named_models_only(  &["--partition-by", "genes", "--annotation", "a.gff", "--model", "jc69", "--shared-params", "frequencies"], "--shared-params applies only with two or more partitions whose model is infer, but --model names a model, whose parameters every partition keeps", &["shared_params", "model"])]
  #[case::one_inferred_entry( &["--partition", "{sites: ['1-10'], model: jc69}", "--partition", "{}", "--shared-params", "frequencies"], "--shared-params applies only with two or more partitions whose model is infer, but the --partition entries have 1 such partition", &["shared_params", "partition"])]
  #[trace]
  fn test_shared_params_argument_errors(#[case] args: &[&str], #[case] message: &str, #[case] settings: &[&str]) {
    let raw = TreetimeOptimizeArgsRaw::try_parse_from([&OPTIMIZE[..], args].concat()).unwrap();
    let expected = (message.to_owned(), settings.iter().map(|setting| (*setting).to_owned()).collect());
    assert_eq!(Some(expected), argument_error(TreetimeOptimizeArgs::try_from(raw)));
  }

  #[test]
  fn test_ancestral_shared_params_with_parsimony_is_an_argument_error() {
    let raw = TreetimeAncestralArgsRaw::try_parse_from(
      [
        &ANCESTRAL[..],
        &["--method-anc", "parsimony", "--shared-params", "none"],
      ]
      .concat(),
    )
    .unwrap();
    let expected = (
      "--shared-params applies only with --method-anc marginal: parsimony reconstruction uses no substitution model"
        .to_owned(),
      vec!["shared_params".to_owned(), "method_anc".to_owned()],
    );
    assert_eq!(Some(expected), argument_error(TreetimeAncestralArgs::try_from(raw)));
  }
}
