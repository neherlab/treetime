#[cfg(test)]
mod tests {
  use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
  use app_output::output_conditions::OutputCondition;
  use app_output::output_plan::OutputSelection;
  use clap::Parser;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::Path;

  #[rustfmt::skip]
  #[rstest]
  #[case::default_set(&["--output-all=out"])]
  #[case::selection(  &["--output-all=out", "--output-selection=nwk,sequences"])]
  #[case::per_file(   &["--output-sequences=seq.fasta"])]
  #[trace]
  fn test_output_requests_auto_mode_plans_sequences_with_a_condition(#[case] outputs: &[&str]) -> Result<(), Report> {
    let args = helpers::timetree(&["--alignment=aln.fasta"], outputs);
    assert_eq!(
      btreemap! { OutputSelection::Sequences => OutputCondition::BranchLengthsReestimated },
      helpers::outputs(args)?.conditional_selections()
    );
    Ok(())
  }

  #[test]
  fn test_output_requests_marginal_mode_plans_sequences_without_a_condition() -> Result<(), Report> {
    let args = helpers::timetree(
      &["--alignment=aln.fasta", "--branch-length-mode=marginal"],
      &["--output-all=out"],
    );
    let resolved = helpers::outputs(args)?;
    assert_eq!(
      (Some(Path::new("out/timetree.sequences.fasta")), btreemap! {}),
      (
        resolved.partition_path(OutputSelection::Sequences, 1),
        resolved.conditional_selections()
      )
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::per_file(
    &["--alignment=aln.fasta", "--branch-length-mode=input"], &["--output-sequences=seq.fasta"],
    ("--output-sequences was requested, but --branch-length-mode input keeps the input branch lengths and reconstructs no sequences",
     &["output_sequences", "branch_length_mode"][..]),
  )]
  #[case::selection(
    &["--alignment=aln.fasta", "--branch-length-mode=input"], &["--output-all=out", "--output-selection=sequences"],
    ("--output-selection lists sequences, but --branch-length-mode input keeps the input branch lengths and reconstructs no sequences",
     &["output_selection", "branch_length_mode"][..]),
  )]
  #[case::no_alignment(
    &["--sequence-length=1000"], &["--output-sequences=seq.fasta"],
    ("--output-sequences was requested, but a run without --alignment reconstructs no sequences",
     &["output_sequences", "alignment"][..]),
  )]
  #[case::table_extension(
    &["--sequence-length=1000"], &["--output-rtt=rtt.ssv"],
    ("The file name 'rtt.ssv' of --output-rtt must end in .tsv (tab-separated) or .csv (comma-separated), optionally \
      followed by a compression suffix such as .gz, or be '-' for standard output",
     &["output_rtt"][..]),
  )]
  #[trace]
  fn test_output_requests_timetree_argument_errors(
    #[case] inputs: &[&str],
    #[case] outputs: &[&str],
    #[case] expected: (&str, &[&str]),
  ) {
    let args = helpers::timetree(inputs, outputs);
    assert_eq!(helpers::owned(expected), helpers::argument_error(helpers::outputs(args)));
  }

  #[test]
  fn test_output_requests_input_mode_leaves_sequences_out_of_all() -> Result<(), Report> {
    let args = helpers::timetree(
      &["--alignment=aln.fasta", "--branch-length-mode=input"],
      &["--output-all=out", "--output-selection=all"],
    );
    assert!(
      !helpers::outputs(args)?
        .paths_by_selection()
        .contains_key(&OutputSelection::Sequences)
    );
    Ok(())
  }

  #[test]
  fn test_output_requests_mugration_confidence_parses_as_a_switch() -> Result<(), Report> {
    let args = TreetimeMugrationArgsRaw::try_parse_from([
      "treetime",
      "--confidence",
      "--tree=t.nwk",
      "--metadata=m.tsv",
      "--attribute=country",
      "--output-all=out",
    ])?;
    assert_eq!(
      Some(Path::new("out/mugration.trait-probabilities.tsv")),
      helpers::outputs(args)?.path(OutputSelection::TraitProbabilities)
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::with_traits(     &["--output-traits=traits.tsv"])]
  #[case::confidence_only( &[])]
  #[trace]
  fn test_output_requests_mugration_confidence_without_output_all_fails(#[case] outputs: &[&str]) {
    let argv = [
      &["treetime", "--confidence", "--tree=t.nwk", "--metadata=m.tsv", "--attribute=country"][..],
      outputs,
    ]
    .concat();
    let args = TreetimeMugrationArgsRaw::try_parse_from(argv).unwrap();
    assert_eq!(
      helpers::owned((
        "--confidence has no effect without --output-all: add --output-all <DIR>, or give the file with \
         --output-trait-probabilities <FILE>",
        &["confidence", "output_trait_probabilities"],
      )),
      helpers::argument_error(helpers::outputs(args))
    );
  }

  mod helpers {
    use crate::command::RawConfig;
    use crate::commands::shared::argument_error::ArgumentError;
    use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
    use app_output::output_plan::ResolvedOutputs;
    use clap::Parser;
    use eyre::Report;
    use itertools::{Itertools, chain};

    pub(super) fn outputs<R: RawConfig>(raw: R) -> Result<ResolvedOutputs, Report> {
      R::wrap(R::Args::try_from(raw)?).resolve_outputs()
    }

    pub(super) fn timetree(inputs: &[&str], outputs: &[&str]) -> TreetimeTimetreeArgsRaw {
      let args = chain!(
        ["treetime", "--tree=t.nwk", "--metadata=m.tsv"],
        inputs.iter().copied(),
        outputs.iter().copied()
      )
      .collect_vec();
      TreetimeTimetreeArgsRaw::try_parse_from(args).expect("timetree arguments parse")
    }

    pub(super) fn argument_error(result: Result<ResolvedOutputs, Report>) -> (String, Vec<String>) {
      let Err(report) = result else {
        panic!("expected an argument error");
      };
      let ArgumentError { message, settings } = report.downcast::<ArgumentError>().expect("expected an argument error");
      (message, settings)
    }

    pub(super) fn owned((message, settings): (&str, &[&str])) -> (String, Vec<String>) {
      (
        message.to_owned(),
        settings.iter().map(|setting| (*setting).to_owned()).collect(),
      )
    }
  }
}
