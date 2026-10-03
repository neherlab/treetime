#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::attach::{complete_alignment_for_leaves, sanitize_to_alphabet};
  use crate::ancestral::partition::AncestralPartition;
  use crate::ancestral::plan::{ReconstructionOptions, ReconstructionPlan, reconstruct_partition};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::Representation;
  use crate::partition::marginal::sample::SampleMode;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
  use crate::test_utils::RecordingSeqSink;
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::{nwk_read_file, nwk_read_str};
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{o, vec_of_owned};

  #[rstest]
  fn test_derived_subs_equal_engine_subs_for_argmax_reconstruction(
    #[values("flu/h3n2/20", "ebola/20", "zika/20")] dataset: &str,
    #[values(
      helpers::fitch(),
      helpers::marginal(Representation::Dense, GtrModelName::Infer),
      helpers::marginal(Representation::Sparse, GtrModelName::Infer)
    )]
    plan: ReconstructionPlan,
    #[values(false, true)] impute: bool,
    #[values(false, true)] include_leaves: bool,
  ) -> Result<(), Report> {
    let (expected, actual) = helpers::engine_and_derived_subs(
      dataset,
      "aln.fasta.xz",
      Alphabet::default(),
      &MutationTrack::Nucleotide,
      plan,
      impute,
      include_leaves,
    )?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rstest]
  fn test_derived_subs_equal_engine_subs_for_amino_acid_cds(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[values(false, true)] impute: bool,
  ) -> Result<(), Report> {
    let (expected, actual) = helpers::engine_and_derived_subs(
      "rsv/a/20",
      "translations/G.fasta.xz",
      Alphabet::new(AlphabetName::Aa)?,
      &MutationTrack::AminoAcid("G".to_owned()),
      helpers::marginal(representation, GtrModelName::Infer),
      impute,
      true,
    )?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rstest]
  #[case::fitch(helpers::fitch())]
  #[case::dense(helpers::marginal(Representation::Dense, GtrModelName::JC69))]
  #[case::sparse(helpers::marginal(Representation::Sparse, GtrModelName::JC69))]
  #[trace]
  fn test_derived_mutations_report_observed_ambiguous_leaf_states(
    #[case] plan: ReconstructionPlan,
    #[values(false, true)] include_leaves: bool,
  ) -> Result<(), Report> {
    let (mutations, sequences) = helpers::ambiguous_leaf_mutations(plan, include_leaves, false)?;
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
      helpers::marginal(representation, GtrModelName::JC69),
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

    type EdgeSubs = BTreeMap<GraphEdgeKey, Vec<Sub>>;

    pub(super) const fn fitch() -> ReconstructionPlan {
      ReconstructionPlan::Fitch
    }

    pub(super) const fn marginal(representation: Representation, model: GtrModelName) -> ReconstructionPlan {
      ReconstructionPlan::Marginal {
        representation,
        model,
        gtr_refinement: None,
      }
    }

    pub(super) fn engine_and_derived_subs(
      dataset: &str,
      alignment: &str,
      alphabet: Alphabet,
      track: &MutationTrack,
      plan: ReconstructionPlan,
      impute: bool,
      include_leaves: bool,
    ) -> Result<(EdgeSubs, EdgeSubs), Report> {
      let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(dataset);
      let parse = nwk_read_file(root.join("tree.nwk"))?;
      let names = parse.names();
      let read_alphabet = match track {
        MutationTrack::Nucleotide => alphabet.clone(),
        MutationTrack::AminoAcid(_) => Alphabet::new(AlphabetName::Aa)?,
      };
      let sequences = read_many_fasta_path(&[root.join(alignment)], &read_alphabet)?
        .into_iter()
        .map(|mut record| {
          record.seq = sanitize_to_alphabet(&record.seq, &alphabet).0;
          AlignmentRecord::from(record)
        })
        .collect();
      let sequences = complete_alignment_for_leaves(&parse.graph, sequences, &alphabet, false, &names, &NoopProgress)?;
      let node_inputs = node_seq_inputs(&parse.graph, &names, sequences);
      let options = ReconstructionOptions::new(include_leaves, impute, SampleMode::Argmax);
      let partition = reconstruct_partition(
        &parse.graph,
        &plan,
        0,
        alphabet,
        &node_inputs,
        &branch_lengths_or_zero(&parse.branch_lengths),
        &options,
        &mut get_random_number_generator(None),
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
            AncestralPartition::Marginal { reconstruction, .. } => reconstruction.edge_subs(graph, edge.key())?,
          };
          Ok((edge.key(), subs.into_iter().sorted_by_key(Sub::pos).collect()))
        })
        .collect::<Result<EdgeSubs, Report>>()?;
      let derived = partition.stream_sequences(graph, track, include_leaves, None)?;
      let actual: EdgeSubs = derived
        .edge_mutations
        .into_iter()
        .map(|(edge_key, mutations)| (edge_key, substitutions(mutations)))
        .collect();
      let alphabet = partition.alphabet();
      let observed_ambiguous = |subs: &[Sub]| -> Vec<usize> {
        subs
          .iter()
          .filter(|sub| !alphabet.is_canonical(sub.qry()))
          .map(Sub::pos)
          .collect()
      };
      let without = |subs: Vec<Sub>, positions: &[usize]| -> Vec<Sub> {
        subs.into_iter().filter(|sub| !positions.contains(&sub.pos())).collect()
      };
      let (expected, actual) = expected
        .into_iter()
        .map(|(edge_key, expected_subs)| {
          let actual_subs = actual[&edge_key].clone();
          let positions = observed_ambiguous(&actual_subs);
          (
            (edge_key, without(expected_subs, &positions)),
            (edge_key, without(actual_subs, &positions)),
          )
        })
        .unzip();
      Ok((expected, actual))
    }

    pub(super) fn ambiguous_leaf_mutations(
      plan: ReconstructionPlan,
      include_leaves: bool,
      impute: bool,
    ) -> Result<(BTreeMap<String, Vec<String>>, BTreeMap<String, String>), Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root;")?;
      let names = parse.names();
      let graph = &parse.graph;
      let sequences = read_many_fasta_str(">A\nAKGANC\n>B\nAGGTAC\n>C\nAGGTAC\n", &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let node_inputs = node_seq_inputs(graph, &names, sequences);
      let partition = reconstruct_partition(
        graph,
        &plan,
        0,
        alphabet,
        &node_inputs,
        &branch_lengths_or_zero(&parse.branch_lengths),
        &ReconstructionOptions::new(include_leaves, impute, SampleMode::Argmax),
        &mut get_random_number_generator(None),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let mut sink = RecordingSeqSink::default();
      let derived = partition.stream_sequences(graph, &MutationTrack::Nucleotide, include_leaves, Some(&mut sink))?;
      let node_sequences: BTreeMap<GraphNodeKey, Seq> =
        sink.items.into_iter().map(|(key, _, seq)| (key, seq)).collect();
      let name = |key: GraphNodeKey| names[&key].clone().expect("all test nodes are named");
      let mut mutations = BTreeMap::new();
      for edge in graph.get_edges() {
        let subs = substitutions(derived.edge_mutations[&edge.key()].clone());
        let mut applied = node_sequences[&edge.source()].clone();
        for sub in &subs {
          applied[sub.pos()] = sub.qry();
        }
        pretty_assertions::assert_eq!(
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
