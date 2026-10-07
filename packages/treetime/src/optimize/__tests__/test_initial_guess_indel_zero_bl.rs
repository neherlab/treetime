#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts};
  use crate::optimize::params::ExistingBranchLengths;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction};
  use crate::seq::indel::InDel;
  use crate::test_utils::dense_reconstruction_mut;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_primitives::seq::Seq;

  const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_auto_preserves_zero_bl_without_indels() -> Result<(), Report> {
    let (graph, reconstruction, mut branch_lengths) = setup_dense(TREE_ZERO_BL)?;

    let total_length = reconstruction.sequence_length();
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
    let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      ExistingBranchLengths::Keep,
      false,
      &mut branch_lengths,
    )?;

    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.key()].unwrap_or(f64::NAN);
      assert!(bl == 0.0, "Without indels, Auto mode should preserve zero BL, got {bl}");
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_auto_overrides_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, mut reconstruction, mut branch_lengths) = setup_dense(TREE_ZERO_BL)?;

    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    {
      let partition = dense_reconstruction_mut(&mut reconstruction);
      partition
        .edges
        .estimates
        .entry(edge_key)
        .or_default()
        .indels
        .push(InDel {
          range: (4, 7),
          seq: Seq::try_from_str("ACG")?,
          kind: crate::seq::indel::InDelKind::Deletion,
        });
    }

    let total_length = reconstruction.sequence_length();
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
    let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      ExistingBranchLengths::Keep,
      false,
      &mut branch_lengths,
    )?;

    let bl = branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap_or(0.0);
    assert!(
      bl > 0.0,
      "Auto mode should override zero BL on indel-bearing edge, got {bl}"
    );

    for edge_ref in graph.get_edges().skip(1) {
      let bl = branch_lengths[&edge_ref.key()].unwrap_or(f64::NAN);
      assert!(bl == 0.0, "Non-indel edge should remain zero, got {bl}");
    }
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn setup_dense(
      newick: &str,
    ) -> Result<(Graph, MarginalReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
        "#}
        .as_bytes(),
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let partition = PartitionMarginalDense::new(alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?;
      let reconstruction =
        MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));

      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

      Ok((graph, reconstruction, branch_lengths))
    }
  }
  use helpers::*;
}
