#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
  use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use helpers::argument_error_of;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{o, vec_of_owned};

  const COALESCENT_SKYLINE: &str = "applies only to --coalescent skyline, but --coalescent is";

  const FITTED_COALESCENT: &str = "applies only to --coalescent opt and --coalescent skyline, but --coalescent is";

  const NEEDS_MARGINAL_MODE: &str =
    "requires re-estimated branch lengths (--branch-length-mode marginal), but --branch-length-mode is input";

  const PARSIMONY_MODEL: &str =
    "applies only with --method-anc marginal: parsimony reconstruction uses no substitution model";

  const ANCESTRAL_MARGINAL: &str = "applies only when a marginal reconstruction runs (--method-anc marginal, or the \
                                    amino-acid reconstruction of --translations), but --method-anc is parsimony";

  #[rustfmt::skip]
  #[rstest]
  #[case::n_skyline_with_opt(              &["--n-skyline", "10", "--coalescent", "opt"],                     Some((format!("--n-skyline {COALESCENT_SKYLINE} opt"),                   ["n_skyline", "coalescent"])))]
  #[case::n_skyline_without_coalescent(    &["--n-skyline", "10"],                                            Some((format!("--n-skyline {COALESCENT_SKYLINE} not set"),               ["n_skyline", "coalescent"])))]
  #[case::stiffness_at_documented_default( &["--skyline-stiffness", "2", "--coalescent", "opt"],              Some((format!("--skyline-stiffness {COALESCENT_SKYLINE} opt"),           ["skyline_stiffness", "coalescent"])))]
  #[case::skyline_settings_with_skyline(   &["--n-skyline", "10", "--skyline-stiffness", "2", "--coalescent", "skyline"], None)]
  #[case::confidence_with_fixed_timescale( &["--coalescent-confidence", "3", "--coalescent", "0.01"],         Some((format!("--coalescent-confidence {FITTED_COALESCENT} 0.01"),       ["coalescent_confidence", "coalescent"])))]
  #[case::confidence_with_opt(             &["--coalescent-confidence", "3", "--coalescent", "opt"],          None)]
  #[case::gen_per_year_without_coalescent( &["--gen-per-year", "50"],                                         Some((o!("--gen-per-year applies only with --coalescent, which is not set"), ["gen_per_year", "coalescent"])))]
  #[case::gen_per_year_with_coalescent(    &["--gen-per-year", "30", "--coalescent", "10y"],                  None)]
  #[case::clock_std_dev_without_intervals( &["--clock-rate", "1e-3", "--clock-std-dev", "1e-4"],              Some((o!("--clock-std-dev applies only with --confidence: only the date intervals use the clock rate uncertainty"), ["clock_std_dev", "confidence"])))]
  #[case::clock_std_dev_with_intervals(    &["--clock-rate", "1e-3", "--clock-std-dev", "1e-4", "--confidence"], None)]
  #[case::partition_rates_one_file(        &["--partition-rates", "equal"],                                   Some((o!("--partition-rates applies only with two or more partitions, but the run has 1 partition: each --alignment file is one partition"), ["partition_rates", "alignment"])))]
  #[case::partition_rates_two_files(       &["--partition-rates", "equal", "--alignment", "b.fasta"],         None)]
  #[case::representation_in_input_mode(    &["--representation", "dense", "--branch-length-mode", "input"],   Some((format!("--representation {NEEDS_MARGINAL_MODE}"),                 ["representation", "branch_length_mode"])))]
  #[case::partition_rates_in_input_mode(   &["--partition-rates", "equal", "--alignment", "b.fasta", "--branch-length-mode", "input"], Some((format!("--partition-rates {NEEDS_MARGINAL_MODE}"), ["partition_rates", "branch_length_mode"])))]
  #[case::threshold_in_input_mode(         &["--ambiguity-threshold", "0.1", "--branch-length-mode", "input"], Some((format!("--ambiguity-threshold {NEEDS_MARGINAL_MODE}"),     ["ambiguity_threshold", "branch_length_mode"])))]
  #[case::threshold_in_marginal_mode(      &["--ambiguity-threshold", "0.1", "--branch-length-mode", "marginal"], None)]
  #[case::representation_in_marginal_mode( &["--representation", "dense", "--branch-length-mode", "marginal"], None)]
  #[case::model_in_input_mode(             &["--model", "jc69", "--model-params", "mu=2", "--branch-length-mode", "input"], None)]
  #[case::sequence_length_with_alignment(  &["--sequence-length", "1000"],                                    Some((o!("--sequence-length applies only without --alignment: the run takes the sequence length from the alignment, and an alignment of variable sites with a separate genome length is not supported"), ["sequence_length", "alignment"])))]
  #[case::documented_defaults(             &[],                                                               None)]
  #[trace]
  fn test_dependent_settings_timetree(#[case] extra: &[&str], #[case] expected: Option<(String, [&str; 2])>) {
    let argv = [&["timetree", "--tree", "t.nwk", "--alignment", "a.fasta"][..], extra].concat();
    assert_eq!(
      expected.map(helpers::owned),
      argument_error_of::<TreetimeTimetreeArgsRaw, TreetimeTimetreeArgs>(&argv)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::cli_zero(    helpers::Source::Cli,    "0",  Some("--n-skyline is 0, but the skyline needs at least 2 grid points"))]
  #[case::cli_one(     helpers::Source::Cli,    "1",  Some("--n-skyline is 1, but the skyline needs at least 2 grid points"))]
  #[case::cli_two(     helpers::Source::Cli,    "2",  None)]
  #[case::config_zero( helpers::Source::Config, "0",  Some("--n-skyline is 0, but the skyline needs at least 2 grid points"))]
  #[case::config_one(  helpers::Source::Config, "1",  Some("--n-skyline is 1, but the skyline needs at least 2 grid points"))]
  #[case::config_two(  helpers::Source::Config, "2",  None)]
  #[trace]
  fn test_dependent_settings_timetree_n_skyline_minimum(
    #[case] source: helpers::Source,
    #[case] n_skyline: &str,
    #[case] expected: Option<&str>,
  ) {
    let raw = helpers::skyline_args(source, n_skyline);
    assert_eq!(
      expected.map(|message| (o!(message), vec_of_owned!["n_skyline"])),
      helpers::conversion_error::<TreetimeTimetreeArgs, _>(raw)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::model_with_parsimony(            &["--method-anc", "parsimony", "--model", "infer"],                Some((format!("--model {PARSIMONY_MODEL}"),                              ["model", "method_anc"])))]
  #[case::model_params_with_parsimony(     &["--method-anc", "parsimony", "--model-params", "mu=1"],          Some((format!("--model-params {PARSIMONY_MODEL}"),                       ["model_params", "method_anc"])))]
  #[case::model_with_parsimony_and_aa(     &["--method-anc", "parsimony", "--model", "jc69", "--translations", "{cds}.fasta"], Some((format!("--model {PARSIMONY_MODEL}"),  ["model", "method_anc"])))]
  #[case::representation_with_parsimony(   &["--method-anc", "parsimony", "--representation", "auto"],        Some((format!("--representation {ANCESTRAL_MARGINAL}"),                  ["representation", "method_anc"])))]
  #[case::threshold_with_parsimony(        &["--method-anc", "parsimony", "--ambiguity-threshold", "0.1"],    Some((format!("--ambiguity-threshold {ANCESTRAL_MARGINAL}"),             ["ambiguity_threshold", "method_anc"])))]
  #[case::threshold_with_parsimony_and_aa( &["--method-anc", "parsimony", "--ambiguity-threshold", "0.1", "--translations", "{cds}.fasta"], None)]
  #[case::representation_with_marginal(    &["--representation", "sparse", "--ambiguity-threshold", "0.1", "--model", "infer"], None)]
  #[case::model_params_without_model(      &["--model-params", "mu=1"],                                       Some((o!("--model-params needs a named --model: infer fits the model from the data and takes no parameters, but got 'mu=1'"), ["model", "model_params"])))]
  #[case::parsimony_with_defaults(         &["--method-anc", "parsimony"],                                    None)]
  #[trace]
  fn test_dependent_settings_ancestral(#[case] extra: &[&str], #[case] expected: Option<(String, [&str; 2])>) {
    let argv = [&["ancestral", "--tree", "t.nwk", "--alignment", "a.fasta"][..], extra].concat();
    assert_eq!(
      expected.map(helpers::owned),
      argument_error_of::<TreetimeAncestralArgsRaw, TreetimeAncestralArgs>(&argv)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::model_with_parsimony(            &["--method-anc", "parsimony", "--model", "jc69"],                 Some((format!("--model {PARSIMONY_MODEL}"),                              ["model", "method_anc"])))]
  #[case::representation_with_parsimony(   &["--method-anc", "parsimony", "--representation", "dense"],       Some((o!("--representation applies only when a marginal reconstruction runs (--method-anc marginal), but --method-anc is parsimony"), ["representation", "method_anc"])))]
  #[case::representation_with_marginal(    &["--representation", "dense", "--model", "jc69"],                None)]
  #[case::parsimony_with_defaults(         &["--method-anc", "parsimony"],                                    None)]
  #[trace]
  fn test_dependent_settings_homoplasy(#[case] extra: &[&str], #[case] expected: Option<(String, [&str; 2])>) {
    let argv = [&["homoplasy", "--tree", "t.nwk", "--alignment", "a.fasta"][..], extra].concat();
    assert_eq!(
      expected.map(helpers::owned),
      argument_error_of::<TreetimeHomoplasyArgsRaw, TreetimeHomoplasyArgs>(&argv)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_merge_shared_mutations_dense(  &["--representation", "dense", "--no-merge-shared-mutations"],  Some((o!("--no-merge-shared-mutations applies only to the sparse representation, but --representation is dense"), ["no_merge_shared_mutations", "representation"])))]
  #[case::no_hoist_reversions_dense(        &["--representation", "dense", "--no-hoist-reversions"],        Some((o!("--no-hoist-reversions applies only to the sparse representation, but --representation is dense"), ["no_hoist_reversions", "representation"])))]
  #[case::no_merge_shared_mutations_sparse( &["--representation", "sparse", "--no-merge-shared-mutations"], None)]
  #[case::no_merge_shared_mutations_auto(   &["--no-merge-shared-mutations", "--no-hoist-reversions"],      None)]
  #[case::partition_rates_one_file(         &["--partition-rates", "estimate"],                             Some((o!("--partition-rates applies only with two or more partitions, but the run has 1 partition: each --alignment file is one partition"), ["partition_rates", "alignment"])))]
  #[case::partition_rates_two_files(        &["--partition-rates", "estimate", "--alignment", "b.fasta"],   None)]
  #[case::documented_defaults(              &["--representation", "dense"],                                 None)]
  #[trace]
  fn test_dependent_settings_optimize(#[case] extra: &[&str], #[case] expected: Option<(String, [&str; 2])>) {
    let argv = [&["optimize", "--tree", "t.nwk", "--alignment", "a.fasta"][..], extra].concat();
    assert_eq!(
      expected.map(helpers::owned),
      argument_error_of::<TreetimeOptimizeArgsRaw, TreetimeOptimizeArgs>(&argv)
    );
  }

  mod helpers {
    use crate::commands::shared::argument_error::ArgumentError;
    use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
    use crate::config::load::{load_config_document, read_settings};
    use crate::config::source::ConfigSource;
    use clap::Parser;
    use eyre::Report;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Source {
      Cli,
      Config,
    }

    pub(super) fn argument_error_of<R: Parser, A: TryFrom<R, Error = Report>>(
      argv: &[&str],
    ) -> Option<(String, Vec<String>)> {
      conversion_error::<A, R>(R::try_parse_from(argv).unwrap())
    }

    pub(super) fn conversion_error<A: TryFrom<R, Error = Report>, R>(raw: R) -> Option<(String, Vec<String>)> {
      A::try_from(raw).err().map(|report| {
        let error = report
          .downcast_ref::<ArgumentError>()
          .unwrap_or_else(|| panic!("expected an argument error, got: {report:?}"));
        (error.message.clone(), error.settings.clone())
      })
    }

    pub(super) fn skyline_args(source: Source, n_skyline: &str) -> TreetimeTimetreeArgsRaw {
      match source {
        Source::Cli => TreetimeTimetreeArgsRaw::try_parse_from([
          "timetree",
          "--tree",
          "t.nwk",
          "--metadata",
          "m.tsv",
          "--alignment",
          "a.fasta",
          "--coalescent",
          "skyline",
          "--n-skyline",
          n_skyline,
        ])
        .unwrap(),
        Source::Config => {
          let text = format!(
            "tree: t.nwk\nmetadata: m.tsv\nalignment: [a.fasta]\ncoalescent: skyline\nn_skyline: {n_skyline}\n"
          );
          let source = ConfigSource::new("timetree.yaml", text.clone());
          let merged = load_config_document::<TreetimeTimetreeArgsRaw>(&source, &text, None).unwrap();
          read_settings::<TreetimeTimetreeArgsRaw>(&merged).unwrap()
        },
      }
    }

    pub(super) fn owned((message, settings): (String, [&str; 2])) -> (String, Vec<String>) {
      (message, settings.map(str::to_owned).to_vec())
    }
  }
}
