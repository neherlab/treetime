#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::cancel::CancelledError;
  use crate::error::OperationError;
  use crate::seq::mutation::{MutationTrack, SequenceMutations, mutation_event_strings, stream_sequence_mutations};
  use crate::seq::sink::{SeqItem, SeqSink};
  use eyre::Report;
  use helpers::{Failure, Stream, error_class, stream};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;
  use treetime_utils::error::report_to_string;
  use treetime_utils::{o, vec_of_owned};

  const TREE: &str = "((A:0.1,B:0.1)X:0.1,C:0.1)root;";

  #[test]
  fn test_mutation_stream_derives_root_and_edge_substitutions() -> Result<(), Report> {
    let sequences = btreemap! {
      "root" => "ACGT",
      "X" => "ACGA",
      "A" => "ACGA",
      "B" => "NCKA",
      "C" => "AC-T",
    };

    let actual = stream(&sequences, None).map_err(OperationError::into_report)?;

    let expected = Stream {
      root_sequence: o!("ACGT"),
      edge_mutations: btreemap! {
        o!("X") => vec_of_owned!["T4A"],
        o!("A") => vec![],
        o!("B") => vec_of_owned!["A1N", "G3K"],
        o!("C") => vec![],
      },
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_mutation_stream_rejects_sequences_of_different_lengths() {
    let sequences = btreemap! {
      "root" => "ACGT",
      "X" => "ACG",
      "A" => "ACG",
      "B" => "ACG",
      "C" => "ACGT",
    };

    let actual = error_class(stream(&sequences, None));

    assert_eq!(
      "InferenceFailed: Parent sequence has length 4, but child sequence has length 3. This is an internal error. Please report it to developers.",
      actual
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::cancelled(Failure::Cancelled,                    "Cancelled")]
  #[case::failed(   Failure::Message("profile unavailable"), "InferenceFailed: profile unavailable")]
  #[trace]
  fn test_mutation_stream_classifies_sequence_errors(#[case] failure: Failure, #[case] expected: &str) {
    let graph = nwk_read_str(TREE).unwrap().graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

    let result = stream_sequence_mutations(
      &graph,
      &alphabet,
      &MutationTrack::Nucleotide,
      true,
      true,
      |_| Err::<Seq, _>(failure.report()),
      |_| vec![],
      None,
    );

    assert_eq!(expected, error_class(result.map(|_| ())));
  }

  #[test]
  fn test_mutation_stream_returns_the_sink_error() {
    let sequences = btreemap! {
      "root" => "ACGT",
      "X" => "ACGT",
      "A" => "ACGT",
      "B" => "ACGT",
      "C" => "ACGT",
    };
    let mut sink = helpers::FailingSink;

    let actual = error_class(stream(&sequences, Some(&mut sink)));

    assert_eq!("SinkFailed: disk full", actual);
  }

  mod helpers {
    use super::*;

    #[derive(Debug, PartialEq)]
    pub(super) struct Stream {
      pub(super) root_sequence: String,
      pub(super) edge_mutations: BTreeMap<String, Vec<String>>,
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Failure {
      Cancelled,
      Message(&'static str),
    }

    impl Failure {
      pub(super) fn report(self) -> Report {
        match self {
          Self::Cancelled => Report::new(CancelledError),
          Self::Message(message) => Report::msg(message),
        }
      }
    }

    pub(super) struct FailingSink;

    impl SeqSink for FailingSink {
      fn emit(&mut self, _item: SeqItem<'_>) -> Result<(), Report> {
        Err(Report::msg("disk full"))
      }
    }

    pub(super) fn stream(
      sequences: &BTreeMap<&str, &str>,
      sink: Option<&mut dyn SeqSink>,
    ) -> Result<Stream, OperationError> {
      let parsed = nwk_read_str(TREE).unwrap();
      let names = parsed.names();
      let graph = parsed.graph;
      let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

      let SequenceMutations {
        root_sequence,
        edge_mutations,
      } = stream_sequence_mutations(
        &graph,
        &alphabet,
        &MutationTrack::Nucleotide,
        true,
        true,
        |key| Seq::try_from_str(sequences[names[&key].as_deref().unwrap()]),
        |_| vec![],
        sink,
      )?;

      let edge_mutations = edge_mutations
        .into_iter()
        .map(|(edge_key, mutations)| {
          let child = graph.get_edge(edge_key).unwrap().target();
          let strings = mutations
            .iter()
            .flat_map(|mutation| mutation_event_strings(&mutation.event).unwrap())
            .collect();
          (names[&child].clone().unwrap(), strings)
        })
        .collect();
      Ok(Stream {
        root_sequence: root_sequence.to_string(),
        edge_mutations,
      })
    }

    pub(super) fn error_class<T>(result: Result<T, OperationError>) -> String {
      match result {
        Ok(_) => o!("Ok"),
        Err(OperationError::Cancelled) => o!("Cancelled"),
        Err(OperationError::InferenceFailed(report)) => format!("InferenceFailed: {}", report_to_string(&report)),
        Err(OperationError::SinkFailed(report)) => format!("SinkFailed: {}", report_to_string(&report)),
        Err(other) => format!("Other: {other}"),
      }
    }
  }
}
