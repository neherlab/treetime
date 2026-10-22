#[cfg(test)]
mod tests {
  use crate::partition::create::Representation;
  use eyre::Report;
  use helpers::{AmbiguitySink, Leaves, reconstructed_sequences};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  const SYMMETRIC_ALIGNMENT: &[(&str, &str)] = &[("A", "AACA"), ("B", "AGTC"), ("C", "AACG"), ("D", "AGTT")];

  #[rustfmt::skip]
  #[rstest]
  #[case::in_order(  "((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)Y:0.1)root;")]
  #[case::permuted(  "((D:0.1,C:0.1)Y:0.1,(B:0.1,A:0.1)X:0.1)root;")]
  #[trace]
  fn test_ambiguous_sites_tied_posteriors_resolve_to_first_state_in_every_leaf_order(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[case] tree: &str,
  ) -> Result<(), Report> {
    let expected = btreemap! {
      "X" => "AACA",
      "Y" => "AACG",
      "root" => "AACA",
    };
    let actual = reconstructed_sequences(representation, tree, SYMMETRIC_ALIGNMENT, Leaves::Omitted, AmbiguitySink::off())?;
    assert_eq!(helpers::owned(&expected), actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::in_order(  "((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)Y:0.1)root;")]
  #[case::permuted(  "((D:0.1,C:0.1)Y:0.1,(B:0.1,A:0.1)X:0.1)root;")]
  #[trace]
  fn test_ambiguous_sites_tied_posteriors_report_iupac_code_of_tied_states(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[case] tree: &str,
  ) -> Result<(), Report> {
    let expected = btreemap! {
      "X" => "ARYM",
      "Y" => "ARYK",
      "root" => "ARYN",
    };
    let actual = reconstructed_sequences(representation, tree, SYMMETRIC_ALIGNMENT, Leaves::Omitted, AmbiguitySink::new(0.0)?)?;
    assert_eq!(helpers::owned(&expected), actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::imputed(  Leaves::Imputed,  btreemap! { "A" => "NA", "B" => "NA", "root" => "NN" })]
  #[case::observed( Leaves::Observed, btreemap! { "A" => "RA", "B" => "YA", "root" => "NN" })]
  #[trace]
  fn test_ambiguous_sites_leaves_report_only_imputed_positions(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[case] leaves: Leaves,
    #[case] expected: BTreeMap<&str, &str>,
  ) -> Result<(), Report> {
    let alignment = &[("A", "RA"), ("B", "YA")];
    let actual = reconstructed_sequences(representation, "(A:0.1,B:0.1)root;", alignment, leaves, AmbiguitySink::new(1.0)?)?;
    assert_eq!(helpers::owned(&expected), actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::off(       None,      btreemap! { "A" => "AA", "B" => "CA", "root" => "AA" })]
  #[case::tie_only(  Some(0.0), btreemap! { "A" => "RA", "B" => "YA", "root" => "NA" })]
  #[trace]
  fn test_ambiguous_sites_imputed_leaf_reports_states_tied_under_its_observation(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[case] threshold: Option<f64>,
    #[case] expected: BTreeMap<&str, &str>,
  ) -> Result<(), Report> {
    let alignment = &[("A", "RA"), ("B", "YA")];
    let sink = threshold.map_or_else(|| Ok(AmbiguitySink::off()), AmbiguitySink::new)?;
    let actual = reconstructed_sequences(representation, "(A:0.1,B:0.1)root;", alignment, Leaves::Imputed, sink)?;
    assert_eq!(helpers::owned(&expected), actual);
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_ambiguous_sites_threshold_one_reports_every_state(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
  ) -> Result<(), Report> {
    let expected = btreemap! {
      "X" => "NNNN",
      "Y" => "NNNN",
      "root" => "NNNN",
    };
    let actual = reconstructed_sequences(
      representation,
      "((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)Y:0.1)root;",
      SYMMETRIC_ALIGNMENT,
      Leaves::Omitted,
      AmbiguitySink::new(1.0)?,
    )?;
    assert_eq!(helpers::owned(&expected), actual);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::partition::create::{Representation, build_marginal_partition};
    use crate::partition::inferred_rate::InferredRate;
    use crate::partition::marginal::posterior_states::AmbiguityThreshold;
    use crate::progress::NoopProgress;
    use crate::seq::sink::{SeqItem, SeqSink, SeqStream, SeqTrack};
    use crate::test_utils::{leaf_seq_inputs, marginal_update_reconstruction};
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) fn reconstructed_sequences(
      representation: Representation,
      tree: &str,
      alignment: &[(&str, &str)],
      leaves: Leaves,
      mut sink: AmbiguitySink,
    ) -> Result<BTreeMap<String, String>, Report> {
      let nwk_parsed = nwk_read(tree.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let fasta = alignment.iter().map(|(name, seq)| format!(">{name}\n{seq}\n")).join("");
      let records: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let reconstruction = build_marginal_partition(
        representation,
        &GtrModel::JC69(JC69Params::default()),
        InferredRate::Unit,
        &graph,
        alphabet,
        leaf_seq_inputs(&graph, &names, records),
        &branch_lengths,
        &NoopProgress,
      )?;
      let (reconstruction, _) = marginal_update_reconstruction(reconstruction, &graph, &branch_lengths)?;
      reconstruction
        .stream_sequences(
          &graph,
          leaves == Leaves::Imputed,
          SeqTrack::Nuc(1),
          leaves != Leaves::Omitted,
          true,
          Some(SeqStream::new(&mut sink, &names)),
        )
        .map_err(|err| err.into_report())?;
      Ok(sink.sequences)
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub(super) enum Leaves {
      Omitted,
      Observed,
      Imputed,
    }

    pub(super) fn owned(expected: &BTreeMap<&str, &str>) -> BTreeMap<String, String> {
      expected
        .iter()
        .map(|(&name, &seq)| (name.to_owned(), seq.to_owned()))
        .collect()
    }

    pub(super) struct AmbiguitySink {
      threshold: Option<AmbiguityThreshold>,
      sequences: BTreeMap<String, String>,
    }

    impl AmbiguitySink {
      pub(super) fn new(threshold: f64) -> Result<Self, Report> {
        Ok(Self {
          threshold: Some(AmbiguityThreshold::new(threshold)?),
          sequences: BTreeMap::new(),
        })
      }

      pub(super) fn off() -> Self {
        Self {
          threshold: None,
          sequences: BTreeMap::new(),
        }
      }
    }

    impl SeqSink for AmbiguitySink {
      fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
        if item.emitted {
          self
            .sequences
            .insert(item.name.to_owned(), item.seq_with_ambiguity().to_string());
        }
        Ok(())
      }

      fn ambiguity(&self) -> Option<AmbiguityThreshold> {
        self.threshold
      }
    }
  }
}
