#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{
    gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts, total_sequence_length,
  };
  use crate::optimize::run_loop::marginal_update_dense;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::indel::InDel;
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_fasta_node_inputs;

  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::seq::Seq;

  /// All-zero branch length tree. Auto mode (overwrite_valid=false) treats
  /// zero BL as "valid" and skips the edge - except when indels are present.
  const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_auto_preserves_zero_bl_without_indels() -> Result<(), Report> {
    let (graph, partitions, mut branch_lengths) = setup_dense(TREE_ZERO_BL)?;

    let total_length = total_sequence_length(&partitions, &[]);
    let indel_counts = gather_edge_indel_counts(&graph, &partitions, &[]);
    let sub_counts = gather_edge_sub_counts(&graph, &partitions, &[])?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &partitions, &[])?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      false,
      false,
      &mut branch_lengths,
    )?;

    // All edges should remain zero: no indels, zero is valid
    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.key()].unwrap_or(f64::NAN);
      assert!(bl == 0.0, "Without indels, Auto mode should preserve zero BL, got {bl}");
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_auto_overrides_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, mut partitions, mut branch_lengths) = setup_dense(TREE_ZERO_BL)?;

    // Inject an indel on the first edge
    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    {
      let partition = &mut partitions[0];
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

    let total_length = total_sequence_length(&partitions, &[]);
    let indel_counts = gather_edge_indel_counts(&graph, &partitions, &[]);
    let sub_counts = gather_edge_sub_counts(&graph, &partitions, &[])?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &partitions, &[])?;
    initial_guess_mixed(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      false,
      false,
      &mut branch_lengths,
    )?;

    // The indel-bearing edge should now have a positive BL
    let bl = branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap_or(0.0);
    assert!(
      bl > 0.0,
      "Auto mode should override zero BL on indel-bearing edge, got {bl}"
    );

    // Non-indel edges should remain zero
    for edge_ref in graph.get_edges().skip(1) {
      let bl = branch_lengths[&edge_ref.key()].unwrap_or(f64::NAN);
      assert!(bl == 0.0, "Non-indel edge should remain zero, got {bl}");
    }
    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn setup_dense(
      newick: &str,
    ) -> Result<(Graph, Vec<DenseReconstruction>, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln = read_many_fasta_str(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
        "#},
        &alphabet,
      )?;
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(&aln)?);
      let node_states = partition.attach_sequences(&graph, &nwk_fasta_node_inputs(&graph, &names, aln))?;
      let partitions = vec![DenseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
        node_states,
      )];

      let (partitions, _) = marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;

      Ok((graph, partitions, branch_lengths))
    }
  }
  use helpers::*;
}
