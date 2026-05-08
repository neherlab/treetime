#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::alignment::get_common_length;
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::optimize_unified::initial_guess_mixed;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::fitch::PartitionFitch;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::traits::PartitionBranchOps;
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGTACGTACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  /// All four sequences share gaps at positions 4-7. Remaining 12 positions
  /// are identical to `gap_free_alignment`.
  fn gappy_alignment_shared() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGT----ACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGT----ACGTACGG
        >D
        ACGT----ACGTACGC
      "#},
      &alphabet,
    )
  }

  /// B alone has gaps at positions 4-7. Other sequences have ACGT there.
  fn gappy_alignment_one_leaf() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  fn setup_sparse(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalSparse>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = PartitionFitch::compress(graph, 0, alphabet, aln)?;
    let partitions = vec![Arc::new(RwLock::new(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?,
    ))];
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

    Ok(partitions)
  }

  fn get_branch_lengths(graph: &GraphAncestral) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect()
  }

  #[test]
  fn test_sparse_effective_length_no_gaps() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse(&graph, &aln)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = partitions[0].read_arc();
      let effective = p.edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_no_gaps() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense(&graph, &aln)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = partitions[0].read_arc();
      let effective = p.edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_shared_gaps() -> Result<(), Report> {
    let aln = gappy_alignment_shared()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse(&graph, &aln)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = partitions[0].read_arc();
      let effective = p.edge_effective_length(&graph, edge_key)?;
      // All nodes share gaps at positions 4-7, so effective = 16 - 4 = 12
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_shared_gaps() -> Result<(), Report> {
    let aln = gappy_alignment_shared()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense(&graph, &aln)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = partitions[0].read_arc();
      let effective = p.edge_effective_length(&graph, edge_key)?;
      // All nodes share gaps at positions 4-7, so effective = 16 - 4 = 12
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_one_leaf_gapped() -> Result<(), Report> {
    let aln = gappy_alignment_one_leaf()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse(&graph, &aln)?;

    let mut found_reduced = false;
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = partitions[0].read_arc();
      let effective = p.edge_effective_length(&graph, edge_key)?;
      // At least one edge (B→AB) should have reduced effective length
      if effective < 16 {
        found_reduced = true;
      }
    }
    assert!(found_reduced);

    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_excludes_gap_positions() -> Result<(), Report> {
    let aln = gappy_alignment_one_leaf()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense(&graph, &aln)?;

    let p = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let subs = p.edge_subs(&graph, edge_key)?;
      // No substitution should involve a gap position (4-7)
      for sub in &subs {
        assert!(sub.pos() < 4 || sub.pos() >= 8);
      }
    }

    Ok(())
  }

  /// With shared gaps at positions 4-7 (which have identical nucleotides in
  /// the gap-free version), the initial guess should produce the same
  /// substitution rate per informative site. Branch lengths should be
  /// proportionally adjusted: subs/12 for gappy vs subs/16 for gap-free.
  #[test]
  fn test_initial_guess_sparse_gap_adjusted_rate() -> Result<(), Report> {
    let aln_clean = gap_free_alignment()?;
    let graph_clean: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions_clean = setup_sparse(&graph_clean, &aln_clean)?;
    initial_guess_mixed(&graph_clean, &partitions_clean, true)?;
    let bl_clean = get_branch_lengths(&graph_clean);

    let aln_gappy = gappy_alignment_shared()?;
    let graph_gappy: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions_gappy = setup_sparse(&graph_gappy, &aln_gappy)?;
    initial_guess_mixed(&graph_gappy, &partitions_gappy, true)?;
    let bl_gappy = get_branch_lengths(&graph_gappy);

    // With 4 shared gap positions out of 16, the effective length is 12.
    // Substitutions at non-gap positions are the same, so the per-site rate
    // is higher by factor 16/12 = 4/3.
    let ratio = 16.0 / 12.0;
    for (clean, gappy) in bl_clean.iter().zip(bl_gappy.iter()) {
      if *clean > 0.0 {
        assert_abs_diff_eq!(gappy / clean, ratio, epsilon = 1e-10);
      } else {
        assert_abs_diff_eq!(*gappy, 0.0, epsilon = 1e-10);
      }
    }

    Ok(())
  }
}
