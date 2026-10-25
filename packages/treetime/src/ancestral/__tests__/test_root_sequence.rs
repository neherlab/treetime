#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::params::{AncestralParams, MethodAncestral};
  use crate::ancestral::pipeline::{AncestralOutput, run};
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::partition::marginal::sample::SampleMode;
  use crate::partition::result::AlignmentSequences;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::SequenceSource;
  use crate::seq::mutation::MutationEvent;
  use crate::test_utils::RecordingSeqSink;
  use crate::test_utils::{ancestral_input_for_tree, leaf_seq_inputs};
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  use std::path::PathBuf;

  use treetime_graph::schedule::GraphSchedule;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::collections::container::get_exactly_one;

  const LONG_BRANCH_TREE: &str = "((A1:0.01,A2:0.01)X:1.0,(C1:0.01,C2:0.01)Y:0.001)root;";
  const LONG_BRANCH_ALIGNMENT: &str = indoc! {r#"
    >A1
    AAAAGGGGAA
    >A2
    AAAAGGGGAA
    >C1
    CAAAGGGGAT
    >C2
    CAAAGGGGAT
  "#};

  #[rstest]
  #[case::sparse(Some(false))]
  #[case::dense(Some(true))]
  fn test_root_sequence_marginal_root_follows_the_short_branch(#[case] dense: Option<bool>) -> Result<(), Report> {
    let (output, _) = helpers::reconstruct(
      LONG_BRANCH_TREE,
      LONG_BRANCH_ALIGNMENT,
      MethodAncestral::Marginal,
      dense,
    )?;
    assert_eq!("CAAAGGGGAT", helpers::sequences(&output).root_sequence.to_string());
    Ok(())
  }

  #[rstest]
  #[case::parsimony(MethodAncestral::Parsimony, None)]
  #[case::sparse(MethodAncestral::Marginal, Some(false))]
  #[case::dense(MethodAncestral::Marginal, Some(true))]
  fn test_root_sequence_equals_augur_root(
    #[case] method: MethodAncestral,
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let (output, emitted) = helpers::reconstruct(LONG_BRANCH_TREE, LONG_BRANCH_ALIGNMENT, method, dense)?;
    assert_eq!(
      emitted[&output.graph.root_key()?],
      helpers::sequences(&output).root_sequence
    );
    Ok(())
  }

  #[rstest]
  #[case::parsimony(MethodAncestral::Parsimony, None)]
  #[case::sparse(MethodAncestral::Marginal, Some(false))]
  #[case::dense(MethodAncestral::Marginal, Some(true))]
  fn test_root_sequence_plus_edge_subs_rebuilds_every_node(
    #[case] method: MethodAncestral,
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let (output, emitted) = helpers::reconstruct(LONG_BRANCH_TREE, LONG_BRANCH_ALIGNMENT, method, dense)?;
    let rebuilt = helpers::rebuild_from_root(&output)?;
    assert_eq!(emitted, rebuilt);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::seq::sink::SeqStream;
    use treetime_graph::key_map::NodeMap;

    pub(super) fn reconstruct(
      nwk: &str,
      fasta: &str,
      method: MethodAncestral,
      dense: Option<bool>,
    ) -> Result<(AncestralOutput, NodeMap<Seq>), Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(nwk.as_bytes())?;
      let names = parse.names.clone();
      let sequences = fasta_read(fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect_vec();
      let length = sequences[0].seq.len();
      let source = SequenceSource {
        path: PathBuf::from("aln.fasta"),
        nodes: leaf_seq_inputs(&parse.graph, &names, sequences),
        alphabet,
        length,
        mask: vec![false; length],
      };
      let input = ancestral_input_for_tree(parse, vec![source], &GtrModel::JC69(JC69Params::default()))?;
      let params = AncestralParams {
        method,
        dense,
        include_leaves: true,
        report_ambiguous: true,
        impute_missing_data: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      };
      let mut sink = RecordingSeqSink::default();
      let output = run(
        &params,
        input,
        Some(SeqStream::new(&mut sink, &names)),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let emitted = sink.items.into_iter().map(|(key, _, seq)| (key, seq)).collect();
      Ok((output, emitted))
    }

    pub(super) fn sequences(output: &AncestralOutput) -> &AlignmentSequences {
      let [result] = output.alignments.as_slice() else {
        panic!("one alignment gives one result");
      };
      result
        .sequences
        .as_ref()
        .expect("ancestral reconstruction returns sequences")
    }

    pub(super) fn rebuild_from_root(output: &AncestralOutput) -> Result<NodeMap<Seq>, Report> {
      let sequences = sequences(output);
      let schedule = GraphSchedule::new(&output.graph)?;
      let mut rebuilt: NodeMap<Seq> = NodeMap::new();
      for &key in schedule.preorder() {
        let parents = schedule.parents(key);
        let seq = if parents.is_empty() {
          sequences.root_sequence.clone()
        } else {
          let (parent_key, edge_key) = get_exactly_one(parents)?;
          let mut seq: Seq = rebuilt[parent_key].clone();
          for mutation in &sequences.edge_mutations[edge_key] {
            if let MutationEvent::Substitution(sub) = &mutation.event {
              seq[sub.pos()] = sub.qry();
            }
          }
          seq
        };
        rebuilt.insert(key, seq);
      }
      Ok(rebuilt)
    }
  }
}
