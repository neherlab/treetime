#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;

  use crate::ancestral::params::AncestralParams;
  use crate::ancestral::params::MethodAncestral;
  use crate::ancestral::pipeline;
  use crate::cancel::NoopCancel;
  use crate::error::OperationError;
  use crate::gtr::get_gtr::GtrModel;
  use crate::partition::marginal::sample::SampleMode;
  use crate::progress::NoopProgress;

  use crate::test_utils::RecordingSeqSink;
  use crate::test_utils::{TestDataset, ancestral_input_for_tree, complete_sequence_source};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;

  #[test]
  fn test_smoke_ancestral_sample_from_profile_root_reproducible() -> Result<(), Report> {
    let seqs_a = helpers::run_root_sampled(42)?;
    let seqs_b = helpers::run_root_sampled(42)?;
    assert_eq!(seqs_a, seqs_b);
    Ok(())
  }

  #[test]
  fn test_sample_from_profile_rejected_for_parsimony() {
    let alphabet = Alphabet::default();
    let parse = nwk_read(TestDataset::Flu20.tree().as_bytes()).unwrap();
    let sequences = TestDataset::Flu20.alignment().unwrap();

    let params = AncestralParams {
      method: MethodAncestral::Parsimony,
      model: GtrModel::Infer,
      dense: None,
      include_leaves: false,
      report_ambiguous: true,
      impute_missing_data: false,
      seed: 0,
      sample_from_profile: SampleMode::Root,
    };
    let names = parse.names.clone();
    let source = complete_sequence_source(&parse.graph, &names, sequences, &alphabet, "aln.fasta").unwrap();
    let input = ancestral_input_for_tree(parse, vec![source]);

    let result = pipeline::run(&params, input, None, &NoopCancel, &NoopProgress, &NoopProgress)
      .map_err(OperationError::into_report);
    assert_error!(
      result,
      "--sample-from-profile=Root requires --method-anc=marginal. Posterior sampling is only defined for marginal reconstruction; Parsimony has no posterior profile to sample."
    );
  }

  #[test]
  fn test_smoke_ancestral_sample_from_profile_all() -> Result<(), Report> {
    let seqs = helpers::run_sampled(SampleMode::All, 7)?;
    assert!(!seqs.is_empty(), "sampled sequences should not be empty");
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::seq::sink::SeqStream;

    pub(super) fn run_root_sampled(seed: u64) -> Result<BTreeMap<String, String>, Report> {
      run_sampled(SampleMode::Root, seed)
    }

    pub(super) fn run_sampled(mode: SampleMode, seed: u64) -> Result<BTreeMap<String, String>, Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(TestDataset::Flu20.tree().as_bytes())?;
      let sequences = TestDataset::Flu20.alignment()?;

      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        model: GtrModel::Infer,
        dense: None,
        include_leaves: false,
        report_ambiguous: true,
        impute_missing_data: false,
        seed,
        sample_from_profile: mode,
      };
      let names = parse.names.clone();
      let source = complete_sequence_source(&parse.graph, &names, sequences, &alphabet, "aln.fasta")?;
      let input = ancestral_input_for_tree(parse, vec![source]);

      let mut sink = RecordingSeqSink::default();
      pipeline::run(
        &params,
        input,
        Some(SeqStream::new(&mut sink, &names)),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let captured = sink
        .items
        .into_iter()
        .filter(|(_, emitted, _)| *emitted)
        .map(|(key, _, seq)| (names[key].to_owned(), seq.to_string()))
        .collect::<BTreeMap<String, String>>();

      Ok(captured)
    }
  }
}
