#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::__tests__::prop_generators::input::{MarginalTestInput, arb_marginal_input_small};
  use crate::ancestral::params::{AncestralParams, MethodAncestral};
  use crate::ancestral::pipeline::{AncestralOutput, run};
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::marginal::sample::SampleMode;
  use crate::partition::result::PartitionSequences;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::{SequenceSource, pair_leaf_sequences};
  use crate::seq::mutation::MutationEvent;
  use crate::test_utils::{RecordingSeqSink, ancestral_input_for_tree};
  use proptest::prelude::*;
  use proptest::test_runner::TestCaseError;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::Seq;

  const RECONSTRUCTIONS: [(MethodAncestral, Option<bool>, SampleMode); 7] = [
    (MethodAncestral::Parsimony, None, SampleMode::Argmax),
    (MethodAncestral::Marginal, Some(true), SampleMode::Argmax),
    (MethodAncestral::Marginal, Some(true), SampleMode::Root),
    (MethodAncestral::Marginal, Some(true), SampleMode::All),
    (MethodAncestral::Marginal, Some(false), SampleMode::Argmax),
    (MethodAncestral::Marginal, Some(false), SampleMode::Root),
    (MethodAncestral::Marginal, Some(false), SampleMode::All),
  ];

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(40))]

    #[test]
    fn test_prop_output_root_and_mutations_rebuild_streamed_sequences(
      input in arb_marginal_input_small(),
      reconstruction in 0..RECONSTRUCTIONS.len(),
      impute in any::<bool>(),
      include_leaves in any::<bool>(),
      seed in any::<u64>(),
    ) {
      let (method, dense, sample_from_profile) = RECONSTRUCTIONS[reconstruction];
      let params = AncestralParams {
        method,
        model: GtrModelName::JC69,
        dense,
        include_leaves,
        report_ambiguous: true,
        impute_missing_data: impute,
        seed,
        sample_from_profile,
      };
      let (sink, output) =
        helpers::run_recorded(&input, &params).map_err(|error| TestCaseError::fail(format!("{error:?}")))?;
      let graph = &output.graph;
      let [result] = output.results.as_slice() else {
        return Err(TestCaseError::fail("one alignment gives one result"));
      };
      let sequences = result.sequences.as_ref().ok_or_else(|| TestCaseError::fail("no sequences"))?;
      let alphabet = Alphabet::default();
      let streamed: BTreeMap<GraphNodeKey, Seq> =
        sink.items.iter().map(|(key, _, seq)| (*key, seq.clone())).collect();
      let root_key = graph.root_key().map_err(|error| TestCaseError::fail(format!("{error:?}")))?;

      prop_assert_eq!(graph.num_nodes(), sink.items.len());
      prop_assert_eq!(&streamed[&root_key], &sequences.root_sequence);
      for edge in graph.get_edges() {
        let parent = &streamed[&edge.source()];
        let child = &streamed[&edge.target()];
        let rebuilt = helpers::apply_substitutions(parent, sequences, edge.key())?;
        for pos in 0..parent.len() {
          if !alphabet.is_gap(parent[pos]) && !alphabet.is_gap(child[pos]) {
            prop_assert_eq!(child[pos], rebuilt[pos], "edge {:?} position {}", edge.key(), pos);
          }
        }
      }
    }
  }

  mod helpers {
    use super::*;
    use crate::seq::sink::SeqStream;
    use eyre::Report;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_primitives::AlignmentRecord;

    pub(super) fn run_recorded(
      input: &MarginalTestInput,
      params: &AncestralParams,
    ) -> Result<(RecordingSeqSink, AncestralOutput), Report> {
      let nwk_parsed = nwk_read(input.newick.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let alphabet = Alphabet::default();
      let alignment: Vec<AlignmentRecord> = input.alignment.clone();
      let sequences = pair_leaf_sequences(&nwk_parsed.graph, &names, alignment).sequences;
      let length = sequences.common_length()?;
      let mask = sequences.mask(length, &alphabet);
      let source = SequenceSource {
        path: PathBuf::from("aln.fasta"),
        alphabet,
        length,
        nodes: sequences.nodes,
        mask,
      };
      let ancestral_input = ancestral_input_for_tree(nwk_parsed, vec![source]);
      let mut sink = RecordingSeqSink::default();
      let output = run(
        params,
        ancestral_input,
        Some(SeqStream::new(&mut sink, &names)),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )
      .map_err(|error| error.into_report())?;
      Ok((sink, output))
    }

    pub(super) fn apply_substitutions(
      parent: &Seq,
      sequences: &PartitionSequences,
      edge_key: GraphEdgeKey,
    ) -> Result<Seq, TestCaseError> {
      let mut rebuilt = parent.clone();
      for mutation in &sequences.edge_mutations[&edge_key] {
        if let MutationEvent::Substitution(sub) = &mutation.event {
          prop_assert_eq!(
            parent[sub.pos()],
            sub.reff(),
            "substitution {} must start from the parent state",
            sub
          );
          rebuilt[sub.pos()] = sub.qry();
        }
      }
      Ok(rebuilt)
    }
  }
}
