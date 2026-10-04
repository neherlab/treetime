#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::{Representation, build_marginal_partition};
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::overlay::SeqOverlay;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::convergence::sequence_changes::{
    AncestralStateSnapshot, capture_ancestral_states, count_sequence_changes,
  };
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};

  const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";

  #[test]
  fn test_count_sequence_changes_counts_differing_and_extra_positions_of_shared_nodes() -> Result<(), Report> {
    let previous: AncestralStateSnapshot = btreemap! {
      GraphNodeKey(0) => SeqOverlay::from(Seq::try_from_str("ACGT")?),
      GraphNodeKey(1) => SeqOverlay::from(Seq::try_from_str("AAAA")?),
      GraphNodeKey(2) => SeqOverlay::from(Seq::try_from_str("CC")?),
    };
    let current: AncestralStateSnapshot = btreemap! {
      GraphNodeKey(0) => SeqOverlay::from(Seq::try_from_str("TCGA")?),
      GraphNodeKey(1) => SeqOverlay::from(Seq::try_from_str("AAAAAA")?),
      GraphNodeKey(3) => SeqOverlay::from(Seq::try_from_str("GG")?),
    };

    assert_eq!(4, count_sequence_changes(&previous, &current));
    Ok(())
  }

  #[test]
  fn test_count_sequence_changes_of_equal_snapshots_is_zero() -> Result<(), Report> {
    let snapshot: AncestralStateSnapshot = btreemap! {
      GraphNodeKey(0) => SeqOverlay::from(Seq::try_from_str("ACGT")?),
    };

    assert_eq!(0, count_sequence_changes(&snapshot, &snapshot.clone()));
    Ok(())
  }

  #[test]
  fn test_capture_ancestral_states_without_reconstruction_is_empty() -> Result<(), Report> {
    let graph = nwk_read(TREE.as_bytes())?.graph;

    assert_eq!(
      BTreeMap::<GraphNodeKey, Seq>::new(),
      helpers::materialize(&capture_ancestral_states(&graph, &BranchModel::Input))
    );
    Ok(())
  }

  #[test]
  fn test_capture_ancestral_states_gives_every_internal_node_a_full_length_sequence() -> Result<(), Report> {
    let (graph, branch_model) = helpers::sparse_model(indoc! {"
      >A
      ACGTACGTAC
      >B
      ACGTACGTAA
      >C
      TCGTACGAAC
      >D
      TCGTACGAAC
    "})?;

    let actual = capture_ancestral_states(&graph, &branch_model);

    let internal: BTreeSet<GraphNodeKey> = graph
      .get_nodes()
      .filter(|node| !node.is_leaf())
      .map(|node| node.key())
      .collect();
    assert_eq!(internal, actual.keys().copied().collect::<BTreeSet<_>>());
    let lengths: BTreeMap<GraphNodeKey, usize> = actual.iter().map(|(key, seq)| (*key, seq.to_seq().len())).collect();
    let expected_lengths: BTreeMap<GraphNodeKey, usize> = internal.iter().map(|key| (*key, 10)).collect();
    assert_eq!(expected_lengths, lengths);
    Ok(())
  }

  #[test]
  fn test_capture_ancestral_states_of_identical_leaves_is_the_shared_sequence() -> Result<(), Report> {
    let (graph, branch_model) = helpers::sparse_model(indoc! {"
      >A
      ACGTACGTAC
      >B
      ACGTACGTAC
      >C
      ACGTACGTAC
      >D
      ACGTACGTAC
    "})?;

    let actual = capture_ancestral_states(&graph, &branch_model);

    let shared = Seq::try_from_str("ACGTACGTAC")?;
    let expected: BTreeMap<GraphNodeKey, Seq> = graph
      .get_nodes()
      .filter(|node| !node.is_leaf())
      .map(|node| (node.key(), shared.clone()))
      .collect();
    assert_eq!(expected, helpers::materialize(&actual));
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn materialize(snapshot: &AncestralStateSnapshot) -> BTreeMap<GraphNodeKey, Seq> {
      snapshot.iter().map(|(key, seq)| (*key, seq.to_seq())).collect()
    }

    pub(super) fn sparse_model(fasta: &str) -> Result<(Graph, BranchModel), Report> {
      let nwk_parsed = nwk_read(TREE.as_bytes())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let reconstruction = build_marginal_partition(
        Representation::Sparse,
        GtrModelName::JC69,
        &graph,
        0,
        alphabet,
        node_seq_inputs(&graph, &names, aln),
        &branch_lengths,
        &NoopProgress,
      )?;
      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths)?;
      Ok((graph, BranchModel::Marginal(reconstruction)))
    }
  }
}
