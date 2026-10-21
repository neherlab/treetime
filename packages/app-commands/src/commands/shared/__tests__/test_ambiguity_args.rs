#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
  use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
  use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
  use crate::yaml::yaml_read_str;
  use app_output::output_plan::OutputSelection;
  use clap::Parser;
  use clap::error::ErrorKind;
  use eyre::Report;
  use helpers::carriers;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::iter::once;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::off(            &[],                                                                          &[])]
  #[case::every_carrier(  &["--ambiguity-threshold=0.125"],                                             &[(OutputSelection::AugurNodeData, 0.125), (OutputSelection::ReconstructedNucFasta, 0.125), (OutputSelection::ReconstructedAaFasta, 0.125)])]
  #[case::chosen_carrier( &["--ambiguity-threshold=0", "--ambiguity-outputs=reconstructed-aa-fasta"],    &[(OutputSelection::ReconstructedAaFasta, 0.0)])]
  #[trace]
  fn test_ambiguity_args_ancestral_carriers(
    #[case] flags: &[&str],
    #[case] expected: &[(OutputSelection, f64)],
  ) -> Result<(), Report> {
    let args = TreetimeAncestralArgsRaw::try_parse_from(once("treetime").chain(flags.iter().copied()))?;
    assert_eq!(expected.to_vec(), carriers(&args.ambiguity.resolve()?));
    Ok(())
  }

  #[test]
  fn test_ambiguity_args_mugration_carries_every_trait_output_by_default() -> Result<(), Report> {
    let args = TreetimeMugrationArgsRaw::try_parse_from(["treetime", "--ambiguity-threshold=0.25"])?;
    let expected = vec![
      (OutputSelection::Auspice, 0.25),
      (OutputSelection::AugurNodeData, 0.25),
      (OutputSelection::TraitsCsv, 0.25),
    ];
    assert_eq!(expected, carriers(&args.ambiguity.resolve()?));
    Ok(())
  }

  #[test]
  fn test_ambiguity_args_timetree_carries_the_reconstructed_fasta() -> Result<(), Report> {
    let args = TreetimeTimetreeArgsRaw::try_parse_from(["treetime", "--ambiguity-threshold=0.125"])?;
    let expected = vec![(OutputSelection::ReconstructedNucFasta, 0.125)];
    assert_eq!(expected, carriers(&args.ambiguity.resolve()?));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::output_of_another_command( &["--ambiguity-threshold=0.1", "--ambiguity-outputs=reconstructed-nuc-fasta"], ErrorKind::InvalidValue)]
  #[case::threshold_above_one(       &["--ambiguity-threshold=1.5"],                                                ErrorKind::ValueValidation)]
  #[case::threshold_not_a_number(    &["--ambiguity-threshold=half"],                                               ErrorKind::ValueValidation)]
  #[trace]
  fn test_ambiguity_args_mugration_rejects_invalid_flags(#[case] flags: &[&str], #[case] expected: ErrorKind) {
    let actual = TreetimeMugrationArgsRaw::try_parse_from(once("treetime").chain(flags.iter().copied()))
      .err()
      .map(|error| error.kind());
    assert_eq!(Some(expected), actual);
  }

  #[test]
  fn test_ambiguity_args_command_line_outputs_without_threshold_is_an_error() -> Result<(), Report> {
    let args = TreetimeMugrationArgsRaw::try_parse_from(["treetime", "--ambiguity-outputs=auspice"])?;
    assert_error!(
      args.ambiguity.resolve(),
      "--ambiguity-outputs requires --ambiguity-threshold"
    );
    Ok(())
  }

  #[test]
  fn test_ambiguity_args_config_outputs_without_threshold_is_an_error() -> Result<(), Report> {
    let args: TreetimeMugrationArgsRaw = yaml_read_str("ambiguity_outputs: [traits-csv]\n")?;
    assert_error!(
      args.ambiguity.resolve(),
      "--ambiguity-outputs requires --ambiguity-threshold"
    );
    Ok(())
  }

  mod helpers {
    use app_output::ambiguity_outputs::AmbiguityOutputs;
    use app_output::output_plan::OutputSelection;
    use strum::IntoEnumIterator;

    pub(super) fn carriers(ambiguity: &AmbiguityOutputs) -> Vec<(OutputSelection, f64)> {
      OutputSelection::iter()
        .filter_map(|output| {
          ambiguity
            .threshold_for(output)
            .map(|threshold| (output, threshold.get()))
        })
        .collect()
    }
  }
}
