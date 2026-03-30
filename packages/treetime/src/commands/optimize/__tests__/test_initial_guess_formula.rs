#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::optimize_unified::initial_guess_mixed;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::traits::PartitionBranchOps;
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// Verify initial_guess_mixed sets branch_length = edge_subs().len() / edge_effective_length()
  /// for sparse partitions. This is the defining formula: the initial branch length
  /// is the fraction of non-gap positions with substitutions.
  #[test]
  fn test_initial_guess_formula_sparse() -> Result<(), Report> {
    let aln = divergent_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse(&graph, &aln)?;

    initial_guess_mixed(&graph, &partitions)?;

    let p = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let sub_count = p.edge_subs(&graph, edge_key)?.len();
      let effective_length = p.edge_effective_length(&graph, edge_key)?;
      let actual_bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);

      let expected_bl = if effective_length > 0 {
        sub_count as f64 / effective_length as f64
      } else {
        0.0
      };

      assert_abs_diff_eq!(expected_bl, actual_bl, epsilon = 1e-15);
    }
    Ok(())
  }

  /// Same formula verification for dense partitions. This test would have caught
  /// the soft Hamming override that was previously shadowing the sub count.
  #[test]
  fn test_initial_guess_formula_dense() -> Result<(), Report> {
    let aln = divergent_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense(&graph, &aln)?;

    initial_guess_mixed(&graph, &partitions)?;

    let p = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let sub_count = p.edge_subs(&graph, edge_key)?.len();
      let effective_length = p.edge_effective_length(&graph, edge_key)?;
      let actual_bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);

      let expected_bl = if effective_length > 0 {
        sub_count as f64 / effective_length as f64
      } else {
        0.0
      };

      assert_abs_diff_eq!(expected_bl, actual_bl, epsilon = 1e-15);
    }
    Ok(())
  }

  fn divergent_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        AAAACCCCGGGGTTTT
        >B
        CCCCGGGGTTTTAAAA
        >C
        GGGGTTTTAAAACCCC
        >D
        TTTTAAAACCCCGGGG
      "#},
      &alphabet,
    )
  }

  fn setup_sparse(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalSparse>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(graph, &partitions, aln)?;
    update_marginal(graph, &partitions)?;

    Ok(partitions)
  }

  fn setup_dense(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalDense>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(graph, &partitions, aln)?;
    update_marginal(graph, &partitions)?;

    Ok(partitions)
  }
}
