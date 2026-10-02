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
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::sync::random::get_random_number_generator;

  #[rstest]
  fn test_derived_subs_equal_engine_subs_for_argmax_reconstruction(
    #[values("flu/h3n2/20", "ebola/20", "zika/20")] dataset: &str,
    #[values(
      helpers::fitch(),
      helpers::marginal(Representation::Dense),
      helpers::marginal(Representation::Sparse)
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
      helpers::marginal(representation),
      impute,
      true,
    )?;
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use super::*;

    type EdgeSubs = BTreeMap<GraphEdgeKey, Vec<Sub>>;

    pub(super) const fn fitch() -> ReconstructionPlan {
      ReconstructionPlan::Fitch
    }

    pub(super) const fn marginal(representation: Representation) -> ReconstructionPlan {
      ReconstructionPlan::Marginal {
        representation,
        model: GtrModelName::Infer,
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
      let actual = derived
        .edge_mutations
        .into_iter()
        .map(|(edge_key, mutations)| (edge_key, substitutions(mutations)))
        .collect();
      Ok((expected, actual))
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
