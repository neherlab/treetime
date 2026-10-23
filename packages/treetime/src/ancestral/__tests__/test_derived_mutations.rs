#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::attach::sanitize_to_alphabet;
  use crate::ancestral::partition::AncestralPartition;
  use crate::ancestral::plan::{ReconstructionMethod, ReconstructionOptions, reconstruct_partition};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::{GtrModel, GtrModelName};
  use crate::partition::create::Representation;
  use crate::partition::inferred_rate::InferredRate;
  use crate::partition::marginal::sample::SampleMode;
  use crate::progress::NoopProgress;
  use crate::seq::mutation::{Mutation, MutationEvent, Sub};
  use crate::seq::sink::SeqTrack;
  use crate::test_utils::RecordingSeqSink;
  use crate::test_utils::complete_leaf_sequences;
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{RSV_A_20_TRANSLATION_G_XZ, TestDataset, fasta_xz_records};
  use eyre::Report;
  use helpers::Plan;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{o, vec_of_owned};

  #[rstest]
  fn test_derived_subs_equal_engine_subs_for_argmax_reconstruction(
    #[values(TestDataset::Flu20, TestDataset::Ebola20, TestDataset::Zika20)] dataset: TestDataset,
    #[values(
      helpers::fitch(),
      helpers::marginal(Representation::Dense, GtrModelName::Infer),
      helpers::marginal(Representation::Sparse, GtrModelName::Infer)
    )]
    plan: Plan,
    #[values(false, true)] impute: bool,
    #[values(false, true)] include_leaves: bool,
  ) -> Result<(), Report> {
    let subs = helpers::engine_and_derived_subs(
      dataset,
      dataset.alignment_xz(),
      Alphabet::default(),
      SeqTrack::Nuc(1),
      &plan,
      impute,
      include_leaves,
    )?;
    assert_eq!(subs.expected, subs.actual);
    Ok(())
  }

  #[rstest]
  fn test_derived_subs_equal_engine_subs_for_amino_acid_cds(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[values(false, true)] impute: bool,
  ) -> Result<(), Report> {
    let subs = helpers::engine_and_derived_subs(
      TestDataset::RsvA20,
      RSV_A_20_TRANSLATION_G_XZ,
      Alphabet::new(AlphabetName::Aa)?,
      SeqTrack::Aa("G"),
      &helpers::marginal(representation, GtrModelName::Infer),
      impute,
      true,
    )?;
    assert_eq!(subs.expected, subs.actual);
    Ok(())
  }

  #[rstest]
  #[case::fitch(helpers::fitch())]
  #[case::dense(helpers::marginal(Representation::Dense, GtrModelName::JC69))]
  #[case::sparse(helpers::marginal(Representation::Sparse, GtrModelName::JC69))]
  #[trace]
  fn test_derived_mutations_report_observed_ambiguous_leaf_states(
    #[case] plan: Plan,
    #[values(false, true)] include_leaves: bool,
  ) -> Result<(), Report> {
    let (mutations, sequences) = helpers::ambiguous_leaf_mutations(&plan, include_leaves, false)?;
    assert_eq!(
      btreemap! {
        o!("A") => vec_of_owned!["G2K", "T4A", "A5N"],
        o!("AB") => vec![],
        o!("B") => vec![],
        o!("C") => vec![],
      },
      mutations
    );
    assert_eq!("AKGANC", sequences["A"]);
    Ok(())
  }

  #[rstest]
  #[case::dense(Representation::Dense)]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_derived_mutations_report_imputed_leaf_states(
    #[case] representation: Representation,
    #[values(false, true)] include_leaves: bool,
  ) -> Result<(), Report> {
    let (mutations, sequences) = helpers::ambiguous_leaf_mutations(
      &helpers::marginal(representation, GtrModelName::JC69),
      include_leaves,
      true,
    )?;
    assert_eq!(
      btreemap! {
        o!("A") => vec_of_owned!["T4A"],
        o!("AB") => vec![],
        o!("B") => vec![],
        o!("C") => vec![],
      },
      mutations
    );
    assert_eq!("AGGAAC", sequences["A"]);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::seq::sink::SeqStream;
    use pretty_assertions::assert_eq;

    type EdgeSubs = BTreeMap<GraphEdgeKey, Vec<Sub>>;

    pub(super) struct EngineAndDerivedSubs {
      pub expected: EdgeSubs,
      pub actual: EdgeSubs,
    }

    pub(super) type Plan = (ReconstructionMethod, GtrModel);

    pub(super) fn fitch() -> Plan {
      (ReconstructionMethod::Fitch, GtrModel::Infer)
    }

    pub(super) fn marginal(representation: Representation, model: GtrModelName) -> Plan {
      (ReconstructionMethod::Marginal(representation), GtrModel::from(model))
    }

    pub(super) fn engine_and_derived_subs(
      dataset: TestDataset,
      alignment_xz: &[u8],
      alphabet: Alphabet,
      track: SeqTrack<'_>,
      plan: &Plan,
      impute: bool,
      include_leaves: bool,
    ) -> Result<EngineAndDerivedSubs, Report> {
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let names = parse.names.clone();
      let read_alphabet = match track {
        SeqTrack::Nuc(_) => alphabet.clone(),
        SeqTrack::Aa(_) => Alphabet::new(AlphabetName::Aa)?,
      };
      let sequences = fasta_xz_records(alignment_xz, &read_alphabet)?
        .into_iter()
        .map(|mut record| {
          record.seq = sanitize_to_alphabet(&record.seq, &alphabet).0;
          AlignmentRecord::from(record)
        })
        .collect();
      let node_inputs = complete_leaf_sequences(&parse.graph, &names, sequences, &alphabet)?.nodes;
      let observed: BTreeMap<GraphNodeKey, Seq> = node_inputs
        .iter()
        .filter_map(|(&key, input)| input.seq.clone().map(|seq| (key, seq)))
        .collect();
      let options = ReconstructionOptions::new(impute, SampleMode::Argmax);
      let partition = reconstruct_partition(
        &parse.graph,
        plan.0,
        &plan.1,
        InferredRate::Unit,
        alphabet,
        node_inputs,
        &branch_lengths_or_zero(&parse.branch_lengths),
        &options,
        &mut get_random_number_generator(0),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let graph = &parse.graph;
      let expected = graph
        .get_edges()
        .map(|edge| {
          let subs = match &partition {
            AncestralPartition::Fitch(fitch) => fitch.edges[&edge.key()].fitch_subs().to_vec(),
            AncestralPartition::Marginal { reconstruction, .. } => {
              reconstruction.view().edge_subs(graph, edge.key())?
            },
          };
          Ok((edge.key(), subs.into_iter().sorted_by_key(Sub::pos).collect()))
        })
        .collect::<Result<EdgeSubs, Report>>()?;
      let mut sink = RecordingSeqSink::default();
      let derived = partition.stream_sequences(
        graph,
        track,
        include_leaves,
        true,
        Some(SeqStream::new(&mut sink, &names)),
      )?;
      let actual: EdgeSubs = derived
        .edge_mutations
        .into_iter()
        .map(|(edge_key, mutations)| (edge_key, substitutions(mutations)))
        .collect();
      let streamed: BTreeMap<GraphNodeKey, Seq> = sink.items.into_iter().map(|(key, _, seq)| (key, seq)).collect();
      let alphabet = partition.alphabet();
      if !impute {
        assert_observed_leaf_states(graph, alphabet, &observed, &streamed, &expected, &actual)?;
      }
      let columns = ambiguous_columns(&observed, alphabet);
      let outside = |subs: &EdgeSubs| -> EdgeSubs {
        subs
          .iter()
          .map(|(&edge_key, subs)| {
            let kept = subs
              .iter()
              .filter(|sub| !columns.contains(&sub.pos()))
              .cloned()
              .collect();
            (edge_key, kept)
          })
          .collect()
      };
      Ok(EngineAndDerivedSubs {
        expected: outside(&expected),
        actual: outside(&actual),
      })
    }

    fn ambiguous_columns(observed: &BTreeMap<GraphNodeKey, Seq>, alphabet: &Alphabet) -> BTreeSet<usize> {
      observed
        .values()
        .flat_map(|seq| {
          seq
            .iter()
            .enumerate()
            .filter(|(_, state)| !alphabet.is_canonical(**state) && !alphabet.is_gap(**state))
            .map(|(pos, _)| pos)
        })
        .collect()
    }

    fn assert_observed_leaf_states(
      graph: &Graph,
      alphabet: &Alphabet,
      observed: &BTreeMap<GraphNodeKey, Seq>,
      streamed: &BTreeMap<GraphNodeKey, Seq>,
      expected: &EdgeSubs,
      actual: &EdgeSubs,
    ) -> Result<(), Report> {
      for edge in graph.get_edges() {
        let Some(leaf) = observed.get(&edge.target()) else {
          continue;
        };
        let parent = &streamed[&edge.source()];
        for (pos, &state) in leaf.iter().enumerate() {
          if alphabet.is_canonical(state) || alphabet.is_gap(state) || alphabet.is_gap(parent[pos]) {
            continue;
          }
          let derived_qry = actual[&edge.key()].iter().find(|sub| sub.pos() == pos).map(Sub::qry);
          let expected_derived = (parent[pos] != state).then_some(state);
          assert_eq!(expected_derived, derived_qry, "edge {:?} position {pos}", edge.key());
          if alphabet.is_canonical(parent[pos]) {
            let engine_state = expected[&edge.key()]
              .iter()
              .find(|sub| sub.pos() == pos)
              .map_or(parent[pos], Sub::qry);
            assert!(
              alphabet.canonical_states(state)?.contains(engine_state),
              "edge {:?} position {pos}: engine state {engine_state} is not a resolution of observed {state}",
              edge.key()
            );
          }
        }
      }
      Ok(())
    }

    pub(super) fn ambiguous_leaf_mutations(
      plan: &Plan,
      include_leaves: bool,
      impute: bool,
    ) -> Result<(BTreeMap<String, Vec<String>>, BTreeMap<String, String>), Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(b"((A:0.1,B:0.1)AB:0.1,C:0.1)root;".as_slice())?;
      let names = parse.names.clone();
      let graph = &parse.graph;
      let sequences = fasta_read(b">A\nAKGANC\n>B\nAGGTAC\n>C\nAGGTAC\n".as_slice(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let node_inputs = leaf_seq_inputs(graph, &names, sequences);
      let partition = reconstruct_partition(
        graph,
        plan.0,
        &plan.1,
        InferredRate::Unit,
        alphabet,
        node_inputs,
        &branch_lengths_or_zero(&parse.branch_lengths),
        &ReconstructionOptions::new(impute, SampleMode::Argmax),
        &mut get_random_number_generator(0),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let mut sink = RecordingSeqSink::default();
      let derived = partition.stream_sequences(
        graph,
        SeqTrack::Nuc(1),
        include_leaves,
        true,
        Some(SeqStream::new(&mut sink, &names)),
      )?;
      let node_sequences: BTreeMap<GraphNodeKey, Seq> =
        sink.items.into_iter().map(|(key, _, seq)| (key, seq)).collect();
      let name = |key: GraphNodeKey| names[key].to_owned();
      let mut mutations = BTreeMap::new();
      for edge in graph.get_edges() {
        let subs = substitutions(derived.edge_mutations[&edge.key()].clone());
        let mut applied = node_sequences[&edge.source()].clone();
        for sub in &subs {
          applied[sub.pos()] = sub.qry();
        }
        assert_eq!(
          node_sequences[&edge.target()].as_str(),
          applied.as_str(),
          "parent sequence plus the mutations of the edge must give the child sequence of {}",
          name(edge.target())
        );
        mutations.insert(name(edge.target()), subs.iter().map(ToString::to_string).collect());
      }
      let sequences = node_sequences
        .iter()
        .map(|(&key, seq)| (name(key), seq.as_str().to_owned()))
        .collect();
      Ok((mutations, sequences))
    }

    fn substitutions(mutations: Vec<Mutation>) -> Vec<Sub> {
      mutations
        .into_iter()
        .filter_map(|mutation| match mutation.event {
          MutationEvent::Substitution(sub) => Some(sub),
          MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
        })
        .collect()
    }
  }
}
