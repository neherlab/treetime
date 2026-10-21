#[cfg(test)]
mod tests {
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_smoke_ancestral_gtr_iterations_sparse() -> Result<(), Report> {
    let initial = helpers::run_flu_20(None, 0)?;
    let refined = helpers::run_flu_20(None, 3)?;

    assert!(
      refined.mu() > 0.0,
      "mu should be positive after GTR iterations, got {}",
      refined.mu()
    );
    assert_ne!(
      format!("{initial:#?}"),
      format!("{refined:#?}"),
      "GTR iterations must replace the parsimony-inferred GTR"
    );

    Ok(())
  }

  #[test]
  fn test_smoke_ancestral_gtr_iterations_dense() -> Result<(), Report> {
    let initial = helpers::run_flu_20(Some(true), 0)?;
    let refined = helpers::run_flu_20(Some(true), 3)?;

    assert!(
      refined.mu() > 0.0,
      "mu should be positive after GTR iterations, got {}",
      refined.mu()
    );
    assert_ne!(
      format!("{initial:#?}"),
      format!("{refined:#?}"),
      "GTR iterations must replace the parsimony-inferred GTR"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::sparse(None      )]
  #[case::dense( Some(true))]
  #[trace]
  fn test_ancestral_one_alignment_with_inferred_model_has_unit_overall_rate(
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let gtr = helpers::run_flu_20(dense, 0)?;

    pretty_assert_ulps_eq!(
      1.0,
      gtr.mu(),
      max_ulps = 0,
      "v0 sets the overall rate of an inferred model to 1: packages/legacy/treetime/treetime/treeanc.py#L1626-L1631"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::sparse(None      )]
  #[case::dense( Some(true))]
  #[trace]
  fn test_ancestral_gtr_iterations_reestimate_the_overall_rate_from_the_unit_start(
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let gtr = helpers::run_flu_20(dense, 3)?;

    assert!(
      (gtr.mu() - 1.0).abs() > 1e-3,
      "the refinement re-estimates mu with W and pi (kb/decisions/ancestral-iterative-gtr-refinement.md), got mu = {}",
      gtr.mu()
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::test_utils::{TestDataset, ancestral_input_for_tree, complete_sequence_source};

    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::gtr::gtr::GTR;
    use crate::partition::marginal::sample::SampleMode;
    use crate::partition::result::PartitionResult;
    use crate::progress::NoopProgress;
    use eyre::{OptionExt, Report};
    use treetime_io::nwk::nwk_read;

    pub(super) fn run_flu_20(dense: Option<bool>, gtr_iterations: usize) -> Result<GTR, Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(TestDataset::Flu20.tree().as_bytes())?;
      let sequences = TestDataset::Flu20.alignment()?;

      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        model: GtrModelName::Infer,
        dense,
        include_leaves: false,
        report_ambiguous: true,
        impute_missing_data: false,
        gtr_iterations,
        site_specific_gtr: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      };
      let names = parse.names.clone();
      let source = complete_sequence_source(&parse.graph, &names, sequences, &alphabet, "aln.fasta")?;
      let input = ancestral_input_for_tree(parse, vec![source]);

      let result = crate::ancestral::pipeline::run(&params, input, None, &NoopCancel, &NoopProgress, &NoopProgress)?;
      let [result] = <[PartitionResult; 1]>::try_from(result.results).expect("one alignment gives one result");
      result.gtr.ok_or_eyre("GTR should be fitted with --model=infer")
    }
  }
}
